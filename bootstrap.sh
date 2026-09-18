#!/usr/bin/env bash
set -euo pipefail
study_checkout=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
study_venv=${STUDY_VENV:-/projects/br1/pzz1/segmentation-sensitivity/software/venv}
module purge
module load GCC/14.3.0 OpenMPI/5.0.8 FEniCS-DOLFINx-Python/0.10.0.post5
python -m venv --system-site-packages "$study_venv"
source "$study_venv/bin/activate"
python -m pip install --no-cache-dir -r "$study_checkout/requirements.txt"
python -m pip install --no-deps --no-build-isolation -e "$study_checkout"
