"""Export a self-contained VTK solution and a report for the supplied network."""
from experiments.synthetic_breaks.config import WORK as SYNTHETIC_WORK
from pathlib import Path
import csv
import hashlib
import html
import json
import shutil
import xml.etree.ElementTree as ET

import h5py
import meshio
import numpy as np
from scipy.spatial import cKDTree

from .run import BASE,load_graph


def read_solution(path,field_name,flux_name):
    with h5py.File(path,'r') as h:
        return (h['/Mesh/mesh/geometry'][:],h['/Mesh/mesh/topology'][:],
                h[f'/Function/{field_name}/0'][:].ravel(),h[f'/Function/{flux_name}/0'][:])


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def verify_xdmf(path):
    for item in ET.parse(path).findall('.//DataItem'):
        file,dataset=item.text.strip().split(':',1)
        with h5py.File(path.parent/file,'r') as h:
            value=h[dataset][:]
        assert value.shape==tuple(int(s) for s in item.get('Dimensions').split()),path
        assert np.isfinite(value).all(),path


def main():
    meta=json.loads((BASE/'input/metadata.json').read_text())
    metrics={}
    for name in ['intact','network_refinement']:
        path=BASE/name
        assert json.loads((path/'complete.json').read_text())['status']=='complete'
        metrics[name]=json.loads((path/'metrics.json').read_text())
        assert json.loads((path/'coverage.json').read_text())['outside_tissue_points']==0
        for kind in ['network','tissue']:verify_xdmf(path/f'{kind}.xdmf')
    main=metrics['intact'];fine=metrics['network_refinement']
    refinement_change=abs(fine['inlet_flux']/main['inlet_flux']-1)
    assert refinement_change<.01,f'Refinement changes inlet flux by {refinement_change:.1%}'
    validation=SYNTHETIC_WORK/'validation'
    old=json.loads((validation/'mpi1/line_n8/intact/metrics.json').read_text())
    regression=json.loads((validation/'import_regression/line_n8/intact/metrics.json').read_text())
    keys=['inlet_flux','wall_exchange','tissue_mean_pressure','tissue_l2_norm']
    regression_error=max(abs(regression[key]/old[key]-1) for key in keys)
    assert regression_error<1e-10
    case=BASE/'intact'
    x3,tets,u3,q3=read_solution(case/'tissue.h5','tissue_pressure','tissue_flux_vector')
    x1,lines,u1,q1=read_solution(case/'network.h5','vessel_pressure','axial_flow_vector')
    graph=load_graph(main['maximum_network_cell_length'])
    graph_xyz=np.array([graph.nodes[i]['pos'] for i in graph])
    distance,point_ids=cKDTree(graph_xyz).query(x1)
    assert distance.max()<1e-12 and len(np.unique(point_ids))==len(graph)
    expected_edges={tuple(sorted(edge)) for edge in graph.edges}
    assert {tuple(sorted(edge)) for edge in point_ids[lines]}==expected_edges
    original_ids=np.array([graph.nodes[int(i)]['original_id'] for i in point_ids],dtype=np.int64)
    inlet=np.asarray(meta['inlet_dimensionless_coordinates'])
    inlet_mask=np.linalg.norm(x1-inlet,axis=1)<1e-10
    assert inlet_mask.sum()==1 and abs(u1[inlet_mask][0]-1)<1e-12
    with np.load(BASE/'input/prepared.npz') as prepared:
        distance,tissue_ids=cKDTree(prepared['tissue']).query(x3)
        assert distance.max()<1e-12 and len(np.unique(tissue_ids))==meta['original_tissue_vertices']
        expected={tuple(sorted(tet)) for tet in prepared['tetrahedra']}
    assert {tuple(sorted(tet)) for tet in tissue_ids[tets]}==expected
    assert np.isclose(main['tissue_volume']*meta['length_scale_m']**3,meta['tissue_volume_m3'],rtol=1e-12)
    cell_data=np.genfromtxt(case/'network_cells.csv',delimiter=',',names=True)
    centers=np.column_stack([cell_data[k] for k in 'xyz'])
    distance,order=cKDTree(centers).query(x1[lines].mean(axis=1))
    assert distance.max()<1e-12 and len(np.unique(order))==len(lines)
    source_edge=cell_data['branch'][order].astype(np.int64)
    directed=[tuple(edge) if graph.has_edge(*edge) else tuple(edge[::-1]) for edge in point_ids[lines]]
    phantom=np.array([graph.edges[edge]['is_phantom'] for edge in directed],dtype=np.int8)
    tangents=np.array([graph_xyz[b]-graph_xyz[a] for a,b in directed])
    tangents/=np.linalg.norm(tangents,axis=1)[:,None]
    assert np.allclose((q1*tangents).sum(axis=1),cell_data['signed_axial_flow'][order],rtol=1e-10,atol=1e-15)
    origin=np.asarray(meta['origin_m']);scale=meta['length_scale_m']
    points_mm=np.vstack([x3,x1])*scale*1000+origin*1000
    n3,n1=len(x3),len(x1);m3,m1=len(tets),len(lines)
    point_data=dict(potential_dimensionless=np.r_[u3,u1],
                    compartment_label=np.r_[np.ones(n3,dtype=np.int8),np.full(n1,2,dtype=np.int8)],
                    is_inlet=np.r_[np.zeros(n3,dtype=np.int8),inlet_mask.astype(np.int8)],
                    source_exodus_node_id=np.r_[tissue_ids+1,np.full(n1,-1,dtype=np.int64)],
                    source_centerline_point_id=np.r_[np.full(n3,-1,dtype=np.int64),original_ids])
    fields=dict(compartment_label=[np.ones(m3,dtype=np.int8),np.full(m1,2,dtype=np.int8)],
                tissue_flux_dimensionless=[q3,np.zeros((m1,3))],
                axial_flow_dimensionless=[np.zeros((m3,3)),q1],
                signed_axial_flow_dimensionless=[np.zeros(m3),cell_data['signed_axial_flow'][order]],
                wall_exchange_per_normalized_length=[np.zeros(m3),cell_data['wall_exchange_per_length'][order]],
                source_centerline_edge_id=[np.full(m3,-1,dtype=np.int64),source_edge],
                radius_mm=[np.zeros(m3),np.full(m1,meta['radius_mm'])],
                is_phantom=[np.zeros(m3,dtype=np.int8),phantom])
    combined=meshio.Mesh(points_mm,[('tetra',tets),('line',lines+n3)],point_data=point_data,cell_data=fields)
    meshio.write(BASE/'coupled_solution.vtu',combined)
    for name,start,stop,block in [('tissue_solution',0,n3,0),('network_solution',n3,n3+n1,1)]:
        meshio.write(BASE/f'{name}.vtu',meshio.Mesh(points_mm[start:stop],
           [(combined.cells[block].type,combined.cells[block].data-start)],
           point_data={key:value[start:stop] for key,value in point_data.items()},
           cell_data={key:[value[block]] for key,value in fields.items()}))
    # Independently read the delivered format and check all fields and geometry.
    for name in ['coupled_solution','network_solution','tissue_solution']:
        loaded=meshio.read(BASE/f'{name}.vtu')
        assert np.isfinite(loaded.points).all()
        for value in loaded.point_data.values():assert np.isfinite(value).all()
        for values in loaded.cell_data.values():
            for value in values:assert np.isfinite(value).all()
        if name=='coupled_solution':
            assert np.allclose(loaded.points,points_mm,rtol=0,atol=1e-12)
            assert np.array_equal(loaded.point_data['potential_dimensionless'],point_data['potential_dimensionless'])
            assert loaded.cells_dict['tetra'].shape==(m3,4) and loaded.cells_dict['line'].shape==(m1,2)
    scalar=[]
    for name,row in metrics.items():scalar.append(dict(run=name,**{key:value for key,value in row.items() if not isinstance(value,(list,dict))}))
    with (BASE/'metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(scalar[0]));writer.writeheader();writer.writerows(scalar)
    checks=dict(status='passed',main_tissue_tetrahedra=m3,main_network_cells=m1,
        original_network_points_preserved=meta['original_network_points'],
        original_tissue_connectivity_preserved=True,network_connectivity_preserved=True,
        coupling_quadrature_inside_tissue=True,inlet_dirichlet_value_verified=True,
        xdmf_hdf5_verified=True,self_contained_vtu_readback_verified=True,
        inlet_flux_refinement_relative_change=refinement_change,
        unchanged_toy_regression_maximum_relative_difference=regression_error)
    (BASE/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    report=f'''# Supplied-network simulation: pt002 inflow

Completed supplement, 11 September 2026. This run solves the intact inflow network from `oilblooddata.zip`, coupled to its supplied tissue mesh. A second solve checks the effect of halving the vessel cell size. Results are dimensionless.

## Completed mesh files

* [coupled_solution.vtu](coupled_solution.vtu): self-contained tissue and vessel mesh with both computed potentials and flux fields. This is the main deliverable and needs no companion file.
* [network_solution.vtu](network_solution.vtu): vessel-only solution, convenient for a Tube filter in ParaView.
* [tissue_solution.vtu](tissue_solution.vtu): full tissue solution.
* [metrics.csv](metrics.csv): the main and refinement measurements.
* [verification.json](verification.json): numerical and artifact checks.

The VTU coordinates are in the original millimetre coordinate frame. `potential_dimensionless` contains tissue u and vessel U. Cell field `compartment_label` is 1 for tissue and 2 for the inflow network; select label 2 with Threshold to expose the embedded network. The axial and tissue flux names end in `dimensionless`; wall exchange is per unit normalized arclength and is positive from vessel to tissue. Signed axial flow follows the source VTP edge orientation. Fields belonging to the other compartment contain zero placeholders. Source VTK point and edge IDs are zero-based; `source_exodus_node_id` preserves the one-based Exodus node IDs. Newly inserted vessel points have source point ID -1. The accompanying native XDMF/HDF5 files retain normalized coordinates.

## Data and inlet selection

The input files are `oilblooddata/thermoembo_run/pt002/pt002_inflow_centerline.vtp` and `oilblooddata/thermoembo_run/pt002/mesh.exo`. The archive README identifies the centerline coordinates as millimetres and the Exodus coordinates as metres. They were converted to one coordinate system before coupling.

The original network has {meta['original_network_points']} vertices, {meta['original_network_edges']} edges, {meta['endpoints']} endpoints, one connected component, and cycle rank {meta['cycle_rank']}. All edges and loops were retained, including {meta['phantom_marked_edges']} edges marked `is_phantom` in the source. Edge subdivision preserves their geometry and connectivity. All stored point and edge radii are 0.1 mm; this supplied value was retained rather than treated as a measured radius distribution.

The inlet is original point **{meta['inlet_original_point_id']}**, at **{meta['inlet_original_coordinates_mm']} mm**. It is the unique degree-one vertex with the highest archived endpoint pressure, {meta['inlet_archived_pressure_mmhg']:.0f} mmHg. This supplies an explicit computational inlet selection. Archived pressures and flows come from a different resistance model and are used only to identify this endpoint; they are not calibration or validation targets for the new PDE solution.

## Model and normalization

The model is the same Laurino-Zunino circular-average 3D-1D weak formulation used in the synthetic pilot:

```text
Integral_Omega K3 grad(u).grad(v)
 + Integral_Lambda K1 A U_prime V_prime
 + Integral_Lambda kappa P (U - T_r u)(V - T_r v) = 0.
```

K3=K1=1, kappa=0.05, A=pi*r^2, and P=2*pi*r. U=1 at the selected inlet; other vessel endpoints have natural zero axial flux. The entire exterior tissue boundary has u=0. Tissue and vessel potentials are continuous P1; circular averages use DG1 with quadrature degree 20. Source terms are zero. The circle checks include every endpoint interpolation location and additional midpoint samples.

Lengths are normalized by the longest extent of the supplied tissue mesh, L={scale*1000:.9f} mm, with origin {meta['origin_m']} metres. Thus r={meta['normalized_radius']:.10g} in the PDE. The imported tissue mesh is preserved: {m3:,} tetrahedra and {n3:,} vertices. Its original volume is {meta['tissue_volume_m3']*1e6:.6f} mL. The main vessel discretization has {m1:,} cells and {n1:,} vertices, with maximum normalized cell length 0.005 (about {scale*1000*.005:.4f} mm). The coupled system has {main['total_dofs']:,} pressure degrees of freedom and was solved on {main['mpi_size']} MPI ranks using PETSc/MUMPS.

These are model potentials and dimensionless fluxes. The coefficients and inlet value have not been fitted to physiological measurements. The explicit use of the supplied radius and organ domain also means the numerical percentages from the earlier synthetic examples should not be transferred to this network.

## Results

'''+table(['Quantity','Computed value'],[
        ['Conservative inlet flux',f"{main['inlet_flux']:.12g}"],
        ['Integrated vessel-to-tissue exchange',f"{main['wall_exchange']:.12g}"],
        ['Tissue exterior outflow',f"{main['tissue_boundary_outflow']:.12g}"],
        ['Mean tissue potential (volume-normalized)',f"{main['tissue_mean_pressure']:.12g}"],
        ['Maximum tissue nodal potential',f'{u3.max():.12g}'],
        ['Vessel nodal potential range',f'{u1.min():.12g} to {u1.max():.12g}'],
        ['Assembly / linear solve time',f"{main['assembly_seconds']:.2f} s / {main['solve_seconds']:.2f} s"]])+f'''

## Validation

'''+table(['Check','Result'],[
        ['Maximum relative flux imbalance',f"{main['conservation_relative_error']:.3e}"],
        ['Relative linear residual',f"{main['linear_relative_residual']:.3e}"],
        ['Relative energy-identity error',f"{main['energy_relative_error']:.3e}"],
        ['Positive inlet value',f"{u1[inlet_mask][0]:.16g}"],
        ['Main coupling containment',f"{main['containment_points']:,} samples; zero outside"],
        ['Halved maximum 1D cell size',f"{fine['network_cells']:,} vessel cells; inlet-flux change {100*refinement_change:.6f}%"],
        ['Existing toy-model regression',f'{regression_error:.3e} maximum relative difference'],
        ['Mesh preservation and VTU read-back','Passed']])+'''

The refinement check concerns the one-dimensional discretization on the supplied tissue mesh. It does not establish three-dimensional mesh convergence. Input checksums, coordinate transformation, source IDs, unmodified extracted inputs, native solutions, and both runs' metrics are included for reproducibility.

## Reproduce

From the study source directory on NOTS:

```bash
source environment.sh
python -m experiments.pt002_inflow.run --prepare
python3 submit.py supplied
# After that job completes:
python -m experiments.pt002_inflow.report
python -m experiments.collect
```

Reference: [Laurino and Zunino (2019), equations (11)-(12)](https://doi.org/10.1051/m2an/2019042).
'''
    (BASE/'REPORT.md').write_text(report)
    # Portable HTML report; geometry is supplied as VTK, without raster rendering.
    import re
    blocks=[];in_table=False;in_code=False
    for line in report.splitlines():
        if line.startswith('```'):
            blocks.append('</pre>' if in_code else '<pre>');in_code=not in_code;continue
        if in_code:blocks.append(html.escape(line)+'\n');continue
        if line.startswith('|'):
            if not in_table:blocks.append('<table>');in_table=True
            cells=[x.strip() for x in line.strip('|').split('|')]
            if all(x=='---' for x in cells):continue
            blocks.append('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in cells)+'</tr>');continue
        if in_table:blocks.append('</table>');in_table=False
        if not line:continue
        escaped=html.escape(line)
        escaped=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',escaped)
        if line.startswith('## '):blocks.append('<h2>'+escaped[3:]+'</h2>')
        elif line.startswith('# '):blocks.append('<h1>'+escaped[2:]+'</h1>')
        else:blocks.append('<p>'+escaped+'</p>')
    style='body{max-width:1050px;margin:35px auto;padding:0 24px;font:16px/1.6 system-ui;color:#172b3a}table{border-collapse:collapse;width:100%}td{padding:8px;border-bottom:1px solid #ddd}tr:first-child{font-weight:600;background:#edf3f6}pre{padding:15px;background:#f4f6f8;overflow:auto}'
    (BASE/'REPORT.html').write_text('<!doctype html><html><meta charset="utf-8"><title>pt002 inflow solution</title><style>'+style+'</style><body>'+''.join(blocks)+'</body></html>')
    print(json.dumps(dict(result=str(BASE),coupled_mesh=str(BASE/'coupled_solution.vtu'),
        mesh_MiB=(BASE/'coupled_solution.vtu').stat().st_size/2**20,
        verification=checks),indent=2))



if __name__=='__main__':main()
