"""Primal P1/P1 discretization of Laurino--Zunino (2019), equations (11)--(12).

The circular average is represented in discontinuous P1 on the network so
different incident branches have independent averaging planes at junctions.
Vessel pressure itself is continuous P1 and satisfies Kirchhoff conservation.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import time

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from scipy.spatial import cKDTree
import ufl
from dolfinx import fem, io, mesh
from fenicsx_ii import Circle, Average, assemble_matrix, create_interpolation_matrix
from .network_mesh import NetworkMesh
from .partition import network_partitioner, tissue_partitioner


@dataclass(frozen=True)
class Parameters:
    radius: float = .025
    kappa: float = .05
    bulk_diffusivity: float = 1.
    axial_diffusivity: float = 1.
    inlet_value: float = 1.
    circle_degree: int = 20

    @property
    def area(self): return np.pi*self.radius**2

    @property
    def exchange(self): return 2*np.pi*self.radius*self.kappa


class Tissue:
    def __init__(self, comm, n, domain=None):
        self.comm = comm
        self.n = n
        self.mesh = domain if domain is not None else mesh.create_unit_cube(comm,n,n,n,partitioner=tissue_partitioner())
        self.V = fem.functionspace(self.mesh,('Lagrange',1))
        self.mesh.topology.create_connectivity(2,3)
        facets = mesh.exterior_facet_indices(self.mesh.topology)
        self.boundary_dofs = fem.locate_dofs_topological(self.V,2,facets)


@dataclass
class Reference:
    tissue: fem.Function
    coordinates: np.ndarray
    pressure: np.ndarray
    metrics: dict


def integral(expr, comm):
    return float(comm.allreduce(fem.assemble_scalar(fem.form(expr)),op=MPI.SUM))


def write_csv(path, names, data):
    with open(path,'w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(names)
        writer.writerows(data)


def solve(tissue, graph, inlet, params, outdir=None, reference=None, volume_output=True):
    comm=tissue.comm
    started=time.perf_counter()
    net=NetworkMesh(graph,N=1,color_strategy='largest_first',comm=comm,
                    cell_partitioner=network_partitioner,build_submeshes=False)
    V=tissue.V
    W=fem.functionspace(net.mesh,('Lagrange',1))
    D=fem.functionspace(net.mesh,('DG',1))
    mixed=ufl.MixedFunctionSpace(V,W)
    u,U=ufl.TrialFunctions(mixed)
    v,w=ufl.TestFunctions(mixed)
    circle=Circle(net.mesh,params.radius,degree=params.circle_degree)
    Tu,Tv=Average(u,circle,D),Average(v,circle,D)
    dx3=ufl.Measure('dx',domain=tissue.mesh)
    dx1=ufl.Measure('dx',domain=net.mesh)
    a=params.bulk_diffusivity*ufl.inner(ufl.grad(u),ufl.grad(v))*dx3
    a+=params.axial_diffusivity*params.area*ufl.inner(ufl.grad(U),ufl.grad(w))*dx1
    a+=params.exchange*(U-Tu)*(w-Tv)*dx1

    nest=assemble_matrix(a)
    indices=[iset.getIndices().copy() for iset in nest.getNestISs()[0]]
    raw=nest.convert('aij')
    raw.assemble()
    A=raw.copy()
    b0=tissue.boundary_dofs
    inlet=np.asarray(inlet).reshape(3,1)
    b1=fem.locate_dofs_geometrical(W,lambda x: np.linalg.norm(x-inlet,axis=0)<1e-10)
    bc0=indices[0][b0[b0<len(indices[0])]]
    bc1=indices[1][b1[b1<len(indices[1])]]
    inlet_count=comm.allreduce(len(bc1))
    if inlet_count != 1:
        raise RuntimeError(f'Expected exactly one inlet DOF, found {inlet_count}')
    lift=A.createVecRight()
    lift.set(0)
    lift.setValues(bc1,np.full(len(bc1),params.inlet_value))
    lift.assemble()
    rhs=A.createVecLeft()
    rhs.set(0)
    bc_all=np.concatenate([bc0,bc1]).astype(PETSc.IntType)
    # PETSc performs the nonzero Dirichlet lifting while eliminating rows AND
    # columns. Preserve the original matrix to recover conservative reactions.
    A.zeroRowsColumns(bc_all,diag=1.,x=lift,b=rhs)
    x=A.createVecRight()
    ksp=PETSc.KSP().create(comm)
    ksp.setOperators(A)
    ksp.setType('preonly')
    ksp.getPC().setType('lu')
    ksp.getPC().setFactorSolverType('mumps')
    ksp.setErrorIfNotConverged(True)
    assembled=time.perf_counter()
    ksp.solve(rhs,x)
    solved=time.perf_counter()
    if ksp.getConvergedReason() <= 0:
        raise RuntimeError(f'Solver did not converge: {ksp.getConvergedReason()}')
    residual=A.createVecLeft()
    A.mult(x,residual)
    residual.axpy(-1,rhs)
    relative_residual=residual.norm()/max(rhs.norm(),1e-30)
    reaction=raw.createVecLeft()
    raw.mult(x,reaction)
    qin=comm.allreduce(float(reaction.getValues(bc1).sum()))
    qout=-comm.allreduce(float(reaction.getValues(bc0).sum()))
    energy=float(x.dot(reaction))
    pressure3=fem.Function(V,name='tissue_pressure')
    pressure1=fem.Function(W,name='vessel_pressure')
    pressure3.x.array[:len(indices[0])]=x.getValues(indices[0])
    pressure1.x.array[:len(indices[1])]=x.getValues(indices[1])
    pressure3.x.scatter_forward()
    pressure1.x.scatter_forward()

    R,_,_=create_interpolation_matrix(V,D,circle,use_petsc=True)
    average=fem.Function(D,name='circular_average_tissue_pressure')
    R.mult(pressure3.x.petsc_vec,average.x.petsc_vec)
    average.x.scatter_forward()
    qexchange=integral(params.exchange*(pressure1-average)*dx1,comm)
    norm3_sq=integral(pressure3**2*dx3,comm)
    norm1_sq=integral(pressure1**2*dx1,comm)
    flow1_sq=integral((params.area*params.axial_diffusivity)**2*
                      ufl.inner(ufl.grad(pressure1),ufl.grad(pressure1))*dx1,comm)
    tissue_volume=integral(ufl.as_ufl(1.)*dx3,comm)
    pressure_integral=integral(pressure3*dx3,comm)
    pressure_mean=pressure_integral/tissue_volume
    conservation=max(abs(qin-qexchange),abs(qin-qout))/max(abs(qin),1e-30)
    energy_error=abs(energy-params.inlet_value*qin)/max(abs(energy),1e-30)
    inlet_error=comm.allreduce(float(np.max(np.abs(x.getValues(bc1)-params.inlet_value),initial=0)),op=MPI.MAX)
    metrics=dict(tissue_n=tissue.n,mpi_size=comm.size,
                 tissue_cells=tissue.mesh.topology.index_map(3).size_global,
                 network_cells=net.mesh.topology.index_map(1).size_global,
                 total_dofs=A.getSize()[0],inlet_flux=qin,wall_exchange=qexchange,
                 tissue_boundary_outflow=qout,tissue_mean_pressure=pressure_mean,
                 tissue_volume=tissue_volume,tissue_integral_pressure=pressure_integral,
                 tissue_l2_norm=np.sqrt(norm3_sq),vessel_l2_norm=np.sqrt(norm1_sq),
                 axial_flow_l2_norm=np.sqrt(flow1_sq),energy=energy,
                 conservation_relative_error=conservation,energy_relative_error=energy_error,
                 linear_relative_residual=relative_residual,inlet_bc_error=inlet_error,
                 assembly_seconds=assembled-started,solve_seconds=solved-assembled,
                 ksp_reason=int(ksp.getConvergedReason()))

    delta3=fem.Function(V,name='tissue_pressure_difference')
    reference1=fem.Function(W,name='intact_pressure_on_retained_vessel')
    delta1=fem.Function(W,name='vessel_pressure_difference')
    coordinates=W.tabulate_dof_coordinates()
    if reference is None:
        reference1.x.array[:]=pressure1.x.array
        metrics.update(relative_inlet_flux_change=0.,relative_tissue_mean_change=0.,
                       tissue_relative_l2_error=0.,vessel_retained_relative_l2_error=0.,
                       axial_flow_retained_relative_l2_error=0.)
    else:
        delta3.x.array[:]=pressure3.x.array-reference.tissue.x.array
        distances,near=cKDTree(reference.coordinates).query(coordinates)
        if np.max(distances,initial=0)>1e-9:
            raise RuntimeError('Perturbed mesh contains nodes absent from intact mesh')
        reference1.x.array[:]=reference.pressure[near]
        reference1.x.scatter_forward()
        delta1.x.array[:]=pressure1.x.array-reference1.x.array
        e3=integral(delta3**2*dx3,comm)
        e1=integral(delta1**2*dx1,comm)
        ref1=integral(reference1**2*dx1,comm)
        eq=integral(ufl.inner(ufl.grad(delta1),ufl.grad(delta1))*dx1,comm)
        refq=integral(ufl.inner(ufl.grad(reference1),ufl.grad(reference1))*dx1,comm)
        metrics.update(relative_inlet_flux_change=qin/reference.metrics['inlet_flux']-1,
                       relative_tissue_mean_change=pressure_mean/reference.metrics['tissue_mean_pressure']-1,
                       tissue_relative_l2_error=np.sqrt(e3)/reference.metrics['tissue_l2_norm'],
                       vessel_retained_relative_l2_error=np.sqrt(e1/max(ref1,1e-30)),
                       axial_flow_retained_relative_l2_error=np.sqrt(eq/max(refq,1e-30)))

    # Gather only the small network for CSV export and exact retained-cell comparisons.
    owned=W.dofmap.index_map.size_local
    point_data=np.column_stack([coordinates[:owned],pressure1.x.array[:owned],
                                reference1.x.array[:owned],delta1.x.array[:owned]])
    points=np.vstack(comm.allgather(point_data))
    flow_space=fem.functionspace(net.mesh,('DG',0,(3,)))
    flow=fem.Function(flow_space,name='axial_flow_vector')
    flow.interpolate(fem.Expression(-params.area*params.axial_diffusivity*ufl.grad(pressure1),
                                     flow_space.element.interpolation_points))
    flow.x.scatter_forward()
    branch_edges=list(graph.edges(data=True))
    centers=np.array([(np.array(graph.nodes[a]['pos'])+np.array(graph.nodes[b]['pos']))/2
                      for a,b,_ in branch_edges])
    edge_tree=cKDTree(centers)
    cell_rows=[]
    local_cells=net.mesh.topology.index_map(1).size_local
    for c in range(local_cells):
        dofs=W.dofmap.cell_dofs(c)
        center=coordinates[dofs].mean(axis=0)
        distance,eidx=edge_tree.query(center)
        if distance>1e-9:
            raise RuntimeError('Could not map a distributed cell to its original edge')
        a,b,attr=branch_edges[eidx]
        pa,pb=np.array(graph.nodes[a]['pos']),np.array(graph.nodes[b]['pos'])
        tangent=(pb-pa)/np.linalg.norm(pb-pa)
        fdof=flow_space.dofmap.cell_dofs(c)[0]
        qvector=flow.x.array[3*fdof:3*fdof+3]
        tissue_avg=average.x.array[D.dofmap.cell_dofs(c)].mean()
        wall=params.exchange*(pressure1.x.array[dofs].mean()-tissue_avg)
        cell_rows.append([attr['original_edge'],attr['branch'],*center,attr['length'],
                          float(qvector@tangent),wall,pressure1.x.array[dofs].mean()])
    cells=comm.gather(cell_rows,root=0)
    if outdir is not None:
        outdir=Path(outdir)
        if comm.rank==0:
            outdir.mkdir(parents=True,exist_ok=True)
            write_csv(outdir/'network_points.csv',
                      ['x','y','z','pressure','intact_pressure','difference'],points)
            write_csv(outdir/'network_cells.csv',
                      ['original_edge','branch','x','y','z','length','signed_axial_flow',
                       'wall_exchange_per_length','mean_pressure'],
                      sorted([r for rows in cells for r in rows],key=lambda r:r[0]))
        comm.barrier()
        with io.XDMFFile(comm,str(outdir/'network.xdmf'),'w') as output:
            output.write_mesh(net.mesh)
            for field in [pressure1,reference1,delta1,flow]: output.write_function(field)
        if volume_output:
            bulk_flow_space=fem.functionspace(tissue.mesh,('DG',0,(3,)))
            bulk_flow=fem.Function(bulk_flow_space,name='tissue_flux_vector')
            bulk_flow.interpolate(fem.Expression(-params.bulk_diffusivity*ufl.grad(pressure3),
                                                  bulk_flow_space.element.interpolation_points))
            with io.XDMFFile(comm,str(outdir/'tissue.xdmf'),'w') as output:
                output.write_mesh(tissue.mesh)
                for field in [pressure3,delta3,bulk_flow]: output.write_function(field)
        if comm.rank==0:
            import meshio
            graph_points=np.array([graph.nodes[n]['pos'] for n in graph.nodes])
            meshio.write(outdir/'geometry.vtu',meshio.Mesh(graph_points,
                         [('line',np.array(list(graph.edges),dtype=np.int64))]))
    metrics['total_seconds']=time.perf_counter()-started
    if not np.isfinite(qin) or qin<=0 or conservation>1e-8 or energy_error>1e-8 or relative_residual>1e-9 or inlet_error>1e-12:
        raise RuntimeError(f'Numerical validation failed: {metrics}')
    # Preserve copies needed by subsequent cases before releasing PETSc resources.
    new_reference=Reference(pressure3,points[:,:3].copy(),points[:,3].copy(),metrics.copy())
    for obj in [ksp,R,residual,reaction,lift,rhs,x,A,raw,nest]:
        obj.destroy()
    return metrics,new_reference
