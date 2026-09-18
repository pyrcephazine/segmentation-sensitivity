#!/usr/bin/env bash
# Source from any directory on NOTS. All project imports come from this checkout.
module purge
module load GCC/14.3.0 OpenMPI/5.0.8 FEniCS-DOLFINx-Python/0.10.0.post5
source "${STUDY_VENV:-/projects/br1/pzz1/segmentation-sensitivity/software/venv}/bin/activate"
study_checkout=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export STUDY_ROOT="$study_checkout"
export STUDY_WORK="${STUDY_WORK:-/scratch/pzz1/segmentation-sensitivity}"
export PYTHONPATH="$study_checkout/src:$study_checkout${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="${MPLCONFIGDIR:-$study_checkout/.work/matplotlib}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$study_checkout/.work/cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
ulimit -c 0
