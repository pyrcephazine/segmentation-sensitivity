import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(120, repeat=False)
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import networkx as nx
import ufl
from dolfinx import fem, mesh
from segmentation_sensitivity.network_mesh import NetworkMesh
from fenicsx_ii import Circle, Average, assemble_matrix
from segmentation_sensitivity.partition import network_partitioner, tissue_partitioner

comm=MPI.COMM_WORLD
G=nx.DiGraph()
for i,xcoord in enumerate(np.linspace(.1,.9,41)):
    G.add_node(i,pos=[float(xcoord),.5,.5])
for i in range(40): G.add_edge(i,i+1)
network=NetworkMesh(G,N=1,color_strategy='largest_first',comm=comm,
                    cell_partitioner=network_partitioner,build_submeshes=False)
m=mesh.create_unit_cube(comm,8,8,8,partitioner=tissue_partitioner())
V=fem.functionspace(m,('Lagrange',1))
W=fem.functionspace(network.mesh,('Lagrange',1))
Wavg=fem.functionspace(network.mesh,('DG',1))
mixed=ufl.MixedFunctionSpace(V,W)
u,U=ufl.TrialFunctions(mixed); v,w=ufl.TestFunctions(mixed)
circle=Circle(network.mesh,.025,degree=20)
Tu=Average(u,circle,Wavg); Tv=Average(v,circle,Wavg)
dx=ufl.Measure('dx',domain=m)
ds=ufl.Measure('dx',domain=network.mesh)
a=ufl.inner(ufl.grad(u),ufl.grad(v))*dx
a+=np.pi*.025**2*ufl.inner(ufl.grad(U),ufl.grad(w))*ds
a+=.05*2*np.pi*.025*(U-Tu)*(w-Tv)*ds
print(comm.rank,'assembling',flush=True)
nest=assemble_matrix(a)
indices=[i.getIndices().copy() for i in nest.getNestISs()[0]]
A=nest.convert('aij'); A.assemble()
raw=A.copy()
boundary=mesh.locate_entities_boundary(m,2,lambda x: np.full(x.shape[1],True))
b0=fem.locate_dofs_topological(V,2,boundary)
b1=fem.locate_dofs_geometrical(W,lambda x: np.linalg.norm(x-np.array([[.1],[.5],[.5]]),axis=0)<1e-10)
bc0=indices[0][b0[b0<len(indices[0])]]
bc1=indices[1][b1[b1<len(indices[1])]]
lift=A.createVecRight(); lift.set(0); lift.setValues(bc1,np.ones(len(bc1))); lift.assemble()
b=A.createVecLeft(); b.set(0)
A.zeroRowsColumns(np.concatenate([bc0,bc1]).astype(PETSc.IntType),diag=1.,x=lift,b=b)
x=A.createVecRight(); ksp=PETSc.KSP().create(comm); ksp.setOperators(A)
ksp.setType('preonly'); ksp.getPC().setType('lu'); ksp.getPC().setFactorSolverType('mumps')
ksp.solve(b,x)
r=raw.createVecLeft(); raw.mult(x,r)
qin=comm.allreduce(float(r.getValues(bc1).sum()))
qout=-comm.allreduce(float(r.getValues(bc0).sum()))
uval=x.getValues(indices[0]); Uval=x.getValues(indices[1])
print(comm.rank,'RESULT',ksp.getConvergedReason(),A.getSize(),'qin',qin,'qout',qout,'network range',Uval.min(),Uval.max(),flush=True)
