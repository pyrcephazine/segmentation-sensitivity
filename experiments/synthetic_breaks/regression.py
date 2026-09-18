"""Run small toy cases and compare with the saved pre-cleanup MPI results.

Use the same arguments as run.py, with STUDY_WORK pointing at a fresh directory.
STUDY_REFERENCE defaults to the historical complete results directory on NOTS.
"""
from pathlib import Path
import json
import os
import sys
from mpi4py import MPI
from .config import WORK, ROOT, RESULTS
from segmentation_sensitivity import network_mesh
from . import run


def main():
    assert Path(network_mesh.__file__).resolve()==ROOT/'src/segmentation_sensitivity/network_mesh.py'
    assert not any(name=='networks_fenicsx' or name.startswith('networks_fenicsx.') for name in sys.modules)
    run.main()
    if MPI.COMM_WORLD.rank!=0:return
    reference=Path(os.environ.get('STUDY_REFERENCE','/home/pzz1/segmentation-sensitivity-results'))/'validation/mpi2'
    checked=[]
    keys=['inlet_flux','wall_exchange','tissue_boundary_outflow','tissue_mean_pressure','tissue_l2_norm','vessel_l2_norm','axial_flow_l2_norm']
    for path in sorted((WORK/'runs').glob('*/*/metrics.json')):
        rel=path.relative_to(WORK/'runs')
        old=json.loads((reference/rel).read_text());new=json.loads(path.read_text())
        used=[key for key in keys if key in old]
        assert len(used)>=5
        errors={key:abs(new[key]-old[key])/max(abs(old[key]),1e-30) for key in used}
        assert max(errors.values())<1e-10,(str(rel),errors)
        assert new['conservation_relative_error']<1e-10
        checked.append(dict(case=str(rel.parent),relative_differences=errors))
    assert len(checked)>=5
    result=dict(status='passed',mpi_ranks=MPI.COMM_WORLD.size,cases=checked,
                local_network_mesh=True,installed_networks_fenicsx_imported=False,
                maximum_relative_difference=max(max(row['relative_differences'].values()) for row in checked))
    RESULTS.mkdir(parents=True,exist_ok=True)
    (RESULTS/'refactor_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
