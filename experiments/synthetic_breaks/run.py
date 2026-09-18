"""Run paired intact/defective cases, recording every seed and removed interval."""
from __future__ import annotations
from .config import WORK
import argparse
import csv
from dataclasses import asdict
import gc
import json
from pathlib import Path
import socket
import sys
import time
import traceback
from mpi4py import MPI
from segmentation_sensitivity.geometry import Geometry
from .cases import cases_for
from segmentation_sensitivity.solver import Parameters, Tissue, solve


def save_json(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def run_suite(args,name):
    comm=MPI.COMM_WORLD
    root=Path(args.output)/f'{name}_n{args.n}'
    geometry={'line':Geometry.line,'y':Geometry.y,'tree':Geometry.tree}[name]()
    all_cases=cases_for(name)
    possible={cut for _,cuts,_ in all_cases for cut in cuts}
    intact=geometry.discretize(args.h,possible)
    cases=all_cases
    if args.cases:
        selected=args.cases.split(',')
        cases=[c for c in cases if c[0]=='intact' or c[0] in selected]
    params=Parameters(inlet_value=args.inlet,circle_degree=args.circle_degree)
    if comm.rank==0:
        root.mkdir(parents=True,exist_ok=True)
        save_json(root/'configuration.json',dict(geometry=name,tissue_n=args.n,
                  maximum_1d_cell_length=args.h,parameters=asdict(params),
                  ranks=comm.size,hostname=socket.gethostname(),
                  branch_points=geometry.points.tolist(),branches=geometry.branches,
                  cases=[dict(name=n,seed=s,breaks=[asdict(b) for b in bs]) for n,bs,s in cases],
                  model=f'Laurino-Zunino 2019 eqs. 11-12; U(inlet)={args.inlet}, u(outer boundary)=0; f=g=0',
                  labels={'inflow':2,'outflow':3},
                  note='Synthetic controlled phase; archive label numbers are not vessel counts.'))
    comm.barrier()
    tissue=Tissue(comm,args.n)
    reference=None
    rows=[]
    for case,cuts,seed in cases:
        graph,stats=geometry.perturb(intact,cuts)
        if comm.rank==0:
            print(f'BEGIN {name} n={args.n} {case}: {len(graph)} network vertices, '
                  f'{stats["n_components"]} components',flush=True)
        output=root/case
        volume=(name!='tree' or case=='intact' or len(cuts)==32 or args.all_volumes)
        metrics,new_reference=solve(tissue,graph,geometry.points[geometry.inlet],params,
                                     output,reference,volume_output=volume)
        row=dict(geometry=name,case=case,seed=seed,**stats,**metrics,
                 saved_tissue_solution=volume)
        rows.append(row)
        if case=='intact': reference=new_reference
        if comm.rank==0:
            save_json(output/'metrics.json',row)
            save_json(root/'metrics.json',rows)
            scalar_rows=[{k:v for k,v in r.items() if not isinstance(v,(list,dict))} for r in rows]
            with open(root/'metrics.csv','w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(scalar_rows[0]))
                writer.writeheader(); writer.writerows(scalar_rows)
            print(f'DONE {case}: Qin={metrics["inlet_flux"]:.9g}, '
                  f'change={100*metrics["relative_inlet_flux_change"]:.3f}%, '
                  f'balance={metrics["conservation_relative_error"]:.2e}, '
                  f'time={metrics["total_seconds"]:.1f}s',flush=True)
        del new_reference
        gc.collect()
        comm.barrier()
    if comm.rank==0:
        save_json(root/'complete.json',dict(completed_cases=len(rows),status='complete',
                  unix_time=time.time(),mpi_size=comm.size))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',choices=['line','y','tree','toys'],default='toys')
    parser.add_argument('--n',type=int,default=24)
    parser.add_argument('--h',type=float,default=.01)
    parser.add_argument('--circle-degree',type=int,default=20)
    parser.add_argument('--inlet',type=float,default=1.)
    parser.add_argument('--cases',help='Comma-separated case names; intact is always included')
    parser.add_argument('--all-volumes',action='store_true')
    parser.add_argument('--output',default=str(WORK/'runs'))
    args=parser.parse_args()
    for name in (['line','y'] if args.geometry=='toys' else [args.geometry]): run_suite(args,name)


if __name__=='__main__':
    try: main()
    except Exception:
        traceback.print_exc()
        sys.stdout.flush(); sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)
