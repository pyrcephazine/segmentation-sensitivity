"""Import and solve one supplied centerline in its matching tissue mesh."""
from pathlib import Path
import argparse
from dataclasses import asdict
import hashlib
import json
import sys
import time
import traceback
import xml.etree.ElementTree as ET
import zipfile

import h5py
import networkx as nx
import numpy as np


from .config import ARCHIVE, WORK
BASE=WORK


def arrays(doc, parent, dtype=float):
    result={}
    for node in doc.findall(f'.//{parent}/DataArray'):
        if node.get('format')!='ascii':raise ValueError('This importer requires ASCII VTP arrays')
        result[node.get('Name')]=np.fromstring(node.text or '',sep=' ',dtype=dtype)
    return result


def prepare():
    inputs=BASE/'input';inputs.mkdir(parents=True,exist_ok=True)
    hashes={}
    with zipfile.ZipFile(ARCHIVE) as archive:
        for name in ['mesh.exo','pt002_inflow_centerline.vtp']:
            member='oilblooddata/thermoembo_run/pt002/'+name
            payload=archive.read(member)
            (inputs/name).write_bytes(payload)
            hashes[member]=hashlib.sha256(payload).hexdigest()
        (inputs/'archive_README.txt').write_bytes(archive.read('oilblooddata/thermoembo_run/README'))
    doc=ET.parse(inputs/'pt002_inflow_centerline.vtp')
    xyz_mm=np.fromstring(doc.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)
    points=arrays(doc,'PointData');data=arrays(doc,'CellData')
    lines=arrays(doc,'Lines',np.int64)
    widths=np.diff(np.r_[0,lines['offsets']])
    if not np.all(widths==2):raise ValueError('Expected two-node centerline cells')
    edges=lines['connectivity'].reshape(-1,2)
    radius_mm=data['radius_mm_edge']
    if not np.allclose(radius_mm,radius_mm[0],rtol=0,atol=1e-14):
        raise ValueError('This single-network run expects the constant radius stored in pt002')
    if radius_mm[0]<=0:raise ValueError('Nonpositive supplied radius')
    graph=nx.Graph();graph.add_nodes_from(range(len(xyz_mm)));graph.add_edges_from(edges)
    if graph.number_of_edges()!=len(edges):raise ValueError('Duplicate supplied edges')
    lengths=np.linalg.norm(xyz_mm[edges[:,1]]-xyz_mm[edges[:,0]],axis=1)
    if lengths.min()<=0:raise ValueError('Zero-length edge')
    if not np.allclose(lengths,data['length_mm'],rtol=1e-10,atol=1e-10):
        raise ValueError('Stored and coordinate-derived edge lengths disagree')
    endpoints=np.array([i for i,degree in graph.degree() if degree==1])
    endpoint_pressures=points['pressure_mmhg'][endpoints]
    maximum=endpoint_pressures.max()
    candidates=endpoints[np.isclose(endpoint_pressures,maximum,rtol=0,atol=1e-8)]
    if len(candidates)!=1:raise ValueError(f'Ambiguous inlet candidates: {candidates}')
    inlet=int(candidates[0])
    with h5py.File(inputs/'mesh.exo') as h:
        if h['connect1'].attrs['elem_type'] not in [b'TETRA',b'TETRA4']:
            raise ValueError('Expected tetrahedra in supplied Exodus mesh')
        tissue_m=h['coord'][:].T.astype(np.float64)
        cells=h['connect1'][:].astype(np.int64)-1
    if cells.shape[1]!=4 or cells.min()!=0 or cells.max()!=len(tissue_m)-1:
        raise ValueError('Unexpected Exodus connectivity')
    origin_m=tissue_m.min(axis=0)
    scale_m=float(np.ptp(tissue_m,axis=0).max())
    normalized_tissue=(tissue_m-origin_m)/scale_m
    normalized_network=(xyz_mm*.001-origin_m)/scale_m
    tet=tissue_m[cells]
    signed_volume=np.linalg.det(np.stack([tet[:,1]-tet[:,0],tet[:,2]-tet[:,0],tet[:,3]-tet[:,0]],axis=2))/6
    if np.any(np.abs(signed_volume)<1e-18):raise ValueError('Degenerate supplied tetrahedron')
    if len(np.unique(tissue_m,axis=0))!=len(tissue_m):raise ValueError('Duplicate supplied tissue vertices')
    if len(np.unique(xyz_mm,axis=0))!=len(xyz_mm):raise ValueError('Duplicate supplied network vertices')
    np.savez_compressed(inputs/'prepared.npz',tissue=normalized_tissue,tetrahedra=cells,
                        network=normalized_network,edges=edges,archived_pressure_mmhg=points['pressure_mmhg'],
                        is_phantom=data['is_phantom'].astype(np.int8),
                        archived_flow_mm3s=data['flow_mm3s'])
    metadata=dict(dataset='pt002',network='inflow',label=2,archive=str(ARCHIVE),input_sha256=hashes,
        original_network_points=len(xyz_mm),original_network_edges=len(edges),
        components=nx.number_connected_components(graph),endpoints=len(endpoints),
        cycle_rank=len(edges)-len(xyz_mm)+nx.number_connected_components(graph),
        maximum_vertex_degree=max(dict(graph.degree()).values()),
        phantom_marked_edges=int(np.count_nonzero(data['is_phantom'])),
        original_tissue_vertices=len(tissue_m),original_tissue_tetrahedra=len(cells),
        tissue_volume_m3=float(np.abs(signed_volume).sum()),
        original_network_length_mm=float(lengths.sum()),radius_mm=float(radius_mm[0]),
        normalized_radius=float(radius_mm[0]*.001/scale_m),
        origin_m=origin_m.tolist(),length_scale_m=scale_m,
        coordinate_transform='x_dimensionless = (x_metres - origin_m) / length_scale_m; VTP millimetres are first multiplied by 0.001',
        inlet_original_point_id=inlet,inlet_original_coordinates_mm=xyz_mm[inlet].tolist(),
        inlet_dimensionless_coordinates=normalized_network[inlet].tolist(),
        inlet_selection='Unique degree-one vertex with highest archived pressure_mmhg among all endpoints',
        inlet_archived_pressure_mmhg=float(maximum),
        archived_pressure_use='Inlet selection only; archived pressures and flows are not boundary data or validation targets for this different model',
        geometry_policy='Retain all original edges, cycles, and phantom-marked edges; only subdivide long cells for the finite element solve')
    (inputs/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2),flush=True)


def load_graph(h):
    with np.load(BASE/'input/prepared.npz') as prepared:
        coords=prepared['network'];edges=prepared['edges'];phantom=prepared['is_phantom']
    graph=nx.DiGraph()
    for i,point in enumerate(coords):graph.add_node(i,pos=point.tolist(),original_id=i)
    cell_id=0
    for original_edge,(a,b) in enumerate(edges):
        a,b=int(a),int(b)
        length=float(np.linalg.norm(coords[b]-coords[a]))
        subdivisions=max(1,int(np.ceil(length/h)))
        ids=[a]
        for k in range(1,subdivisions):
            index=len(graph)
            graph.add_node(index,pos=((1-k/subdivisions)*coords[a]+k/subdivisions*coords[b]).tolist(),original_id=-1)
            ids.append(index)
        ids.append(b)
        for u,v in zip(ids[:-1],ids[1:]):
            graph.add_edge(u,v,original_edge=cell_id,branch=original_edge,
                           length=length/subdivisions,is_phantom=int(phantom[original_edge]))
            cell_id+=1
    return graph


def imported_tissue(comm):
    import basix.ufl
    from dolfinx import mesh
    import ufl
    from segmentation_sensitivity.partition import tissue_partitioner
    from segmentation_sensitivity.solver import Tissue
    if comm.rank==0:
        with np.load(BASE/'input/prepared.npz') as prepared:
            xyz=prepared['tissue'];cells=prepared['tetrahedra']
    else:
        xyz=np.empty((0,3),dtype=np.float64);cells=np.empty((0,4),dtype=np.int64)
    domain=ufl.Mesh(basix.ufl.element('Lagrange','tetrahedron',1,shape=(3,)))
    tissue_mesh=mesh.create_mesh(comm,cells=cells,e=domain,x=xyz,partitioner=tissue_partitioner())
    return Tissue(comm,None,domain=tissue_mesh)


def coverage(tissue,graph,radius,circle_degree,path):
    from dolfinx import geometry
    from fenicsx_ii import Circle
    from segmentation_sensitivity.network_mesh import NetworkMesh
    from mpi4py import MPI
    from segmentation_sensitivity.partition import network_partitioner
    comm=tissue.comm
    net=NetworkMesh(graph,N=1,color_strategy='largest_first',comm=comm,
                    cell_partitioner=network_partitioner,build_submeshes=False)
    circle=Circle(net.mesh,radius,degree=circle_degree)
    n=net.mesh.topology.index_map(1).size_local
    cells=np.arange(n,dtype=np.int32)
    # Endpoints are exactly the DG1 circular-average interpolation locations;
    # midpoint circles add an independent interior containment check.
    quadrature=circle.compute_quadrature(cells,np.array([[0.],[.5],[1.]]))
    points=quadrature.points.reshape(-1,3)
    owners=geometry.determine_point_ownership(tissue.mesh,points,padding=1e-10)
    bad=np.asarray(owners.src_owner)<0
    outside=comm.allreduce(int(bad.sum()),op=MPI.SUM)
    total=comm.allreduce(len(points),op=MPI.SUM)
    samples=comm.gather(points[bad][:5].tolist(),root=0)
    report=dict(total_quadrature_points=total,outside_tissue_points=outside,
                containment_passed=outside==0,sample_outside_points=samples if comm.rank==0 else [])
    if comm.rank==0:
        print('COVERAGE',json.dumps(report),flush=True)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(report,indent=2)+'\n')
    if outside:raise RuntimeError(f'{outside} circular-average points lie outside the supplied tissue')
    return report


def run(args):
    from mpi4py import MPI
    from segmentation_sensitivity.solver import Parameters,solve
    comm=MPI.COMM_WORLD
    metadata=json.loads((BASE/'input/metadata.json').read_text())
    params=Parameters(radius=metadata['normalized_radius'],circle_degree=args.circle_degree)
    graph=load_graph(args.h)
    tissue=imported_tissue(comm)
    out=BASE/args.name
    checked=coverage(tissue,graph,params.radius,params.circle_degree,out/'coverage.json')
    if args.check_only:return
    if comm.rank==0:
        out.mkdir(parents=True,exist_ok=True)
        print(f'BEGIN pt002 inflow: {len(graph)} vertices, {graph.number_of_edges()} cells, {comm.size} MPI ranks',flush=True)
    metrics,reference=solve(tissue,graph,metadata['inlet_dimensionless_coordinates'],params,out,volume_output=True)
    metrics.update(dataset='pt002',geometry='supplied_inflow',case='intact',
                   maximum_network_cell_length=args.h,parameters=asdict(params),
                   topology_components=metadata['components'],topology_cycle_rank=metadata['cycle_rank'],
                   original_input_points=metadata['original_network_points'],
                   original_input_edges=metadata['original_network_edges'],
                   containment_points=checked['total_quadrature_points'])
    if comm.rank==0:
        (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
        (out/'complete.json').write_text(json.dumps(dict(status='complete',unix_time=time.time(),mpi_size=comm.size))+'\n')
        print('DONE',json.dumps(metrics),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--h',type=float,default=.005)
    parser.add_argument('--circle-degree',type=int,default=20)
    parser.add_argument('--name',default='intact')
    args=parser.parse_args()
    if args.prepare:prepare()
    else:run(args)


if __name__=='__main__':
    try:main()
    except Exception:
        traceback.print_exc();sys.stdout.flush();sys.stderr.flush()
        if '--prepare' not in sys.argv:
            from mpi4py import MPI
            MPI.COMM_WORLD.Abort(1)
        raise
