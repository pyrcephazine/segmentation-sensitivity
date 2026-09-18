"""Independent one-dimensional analytic limit of the coupled model."""
from .config import WORK
from dataclasses import asdict
import json
from pathlib import Path
import sys
import traceback
from mpi4py import MPI
import numpy as np
from segmentation_sensitivity.geometry import Geometry
from segmentation_sensitivity.solver import Parameters, Tissue, solve


def main():
    comm=MPI.COMM_WORLD
    root=WORK/'validation/analytic'
    if comm.rank==0: root.mkdir(parents=True,exist_ok=True)
    geometry=Geometry.line()
    tissue=Tissue(comm,8)
    # As bulk diffusivity -> infinity, tissue pressure -> 0. The limiting
    # vessel solution with U(0)=1 and U'(L)=0 has an elementary cosh profile.
    params=Parameters(bulk_diffusivity=1e8)
    length=.8
    ell=np.sqrt(params.axial_diffusivity*params.area/params.exchange)
    exact_q=params.axial_diffusivity*params.area/ell*np.tanh(length/ell)
    rows=[]
    for h in [.04,.02,.01,.005]:
        intact=geometry.discretize(h)
        graph,_=geometry.perturb(intact,[])
        metrics,ref=solve(tissue,graph,geometry.points[0],params)
        s=ref.coordinates[:,0]-.1
        exact=np.cosh((length-s)/ell)/np.cosh(length/ell)
        row=dict(h=h,computed_flux=metrics['inlet_flux'],exact_flux=exact_q,
                 relative_flux_error=abs(metrics['inlet_flux']/exact_q-1),
                 maximum_nodal_pressure_error=float(np.max(np.abs(ref.pressure-exact))),
                 conservation_error=metrics['conservation_relative_error'])
        rows.append(row)
        if comm.rank==0: print('ANALYTIC',row,flush=True)
    errors=np.array([r['relative_flux_error'] for r in rows])
    orders=np.log2(errors[:-1]/errors[1:])
    if not np.all(orders>1.8) or errors[-1]>1e-4:
        raise RuntimeError(f'Analytic convergence test failed: errors={errors}, orders={orders}')
    if comm.rank==0:
        (root/'results.json').write_text(json.dumps(dict(parameters=asdict(params),
             formula='U(s)=cosh((L-s)/ell)/cosh(L/ell); Qin=A/ell*tanh(L/ell)',
             samples=rows,observed_flux_orders=orders.tolist(),status='passed'),indent=2)+'\n')


if __name__=='__main__':
    try: main()
    except Exception:
        traceback.print_exc(); sys.stdout.flush(); sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)
