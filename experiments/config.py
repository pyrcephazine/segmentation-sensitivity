"""Common repository and input paths; no simulation state in the library."""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = Path(os.environ.get('STUDY_DATA_ZIP', '/projects/br1/pzz1/segmentation-sensitivity/oilblooddata.zip'))


def experiment_paths(name):
    results = ROOT / 'experiments' / name / 'results'
    # Published outputs stay with the experiment; large scratch work is optional.
    work = Path(os.environ['STUDY_WORK']) / name if os.environ.get('STUDY_WORK') else results / 'raw'
    return results, work
