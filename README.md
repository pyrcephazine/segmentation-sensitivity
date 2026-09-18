# Segmentation sensitivity

Laurino–Zunino (2019) circular-average 3D–1D diffusion with a positive vessel
inlet. The reusable solver and mesh objects are in `src/segmentation_sensitivity/`.
Experiment scripts and published outputs are in `experiments/<name>/results/`
(outputs) and their parent experiment directory (scripts).

| Experiment | Results |
| --- | --- |
| `synthetic_breaks`: 42 line, Y, and tree cases | [Report](experiments/synthetic_breaks/results/REPORT.html) |
| `pt002_inflow`: supplied network and vessel-resolution check | [Report](experiments/pt002_inflow/results/REPORT.html) · [Mesh](experiments/pt002_inflow/results/coupled_solution.vtu) |

## Submit jobs

One command submits the complete workflow, including reports after all required
simulation and validation jobs succeed:

```bash
python3 submit.py all
```

No manual environment activation or ZIP preparation is needed for submission.
Each job loads `environment.sh`; the supplied-network job prepares its inputs.
Preview the scripts and resource requests without submitting anything:

```bash
python3 submit.py all --dry-run
```

Individual jobs use the same interface:

| Job | Purpose | MPI tasks | Memory | Time limit |
| --- | --- | --- | --- | --- |
| `synthetic` | Main line/Y/tree study and eight-rank check | 8 | 32G | 2h |
| `validate` | Analytic and discretization checks | 4 | 24G | 1h |
| `mpi-check` | Serial/MPI and inlet-amplitude checks | 2 | 16G | 30m |
| `mesh-check` | Large-tree mesh refinement | 4 | 32G | 1h |
| `supplied` | Prepare and solve pt002 at two vessel resolutions | 4 | 16G | 45m |
| `reports` | Generate and publish both experiments' reports | 1 | 8G | 20m |
| `check` | Five-case comparison with historical results | 2 | 8G | 15m |

```bash
python3 submit.py supplied
python3 submit.py reports --after 12345 12346
```

`--after` uses Slurm's `afterok` dependency. The `reports` job needs both completed
experiments and validation data; run it alone only when those inputs exist.
`all` submits the five main jobs independently, then makes reports depend on all
five. It stops on a submission error and prints IDs for jobs already submitted.
Logs go to `.work/logs/<job>-<job-id>.log`. Resource profiles and scientific
commands are defined once in `experiments/jobs.py`; no per-job Slurm files are
needed. Account, partition, and QoS can be overridden with CLI options.

## Environment

The existing NOTS environment is ready to use. Only to recreate it:

```bash
bash bootstrap.sh
```

For interactive commands, use `source environment.sh`. The native FEniCSx 0.10,
PETSc/MUMPS, MPI, and fenicsx_ii 0.4 packages are external and unmodified. Bootstrap
installs the package in editable mode and additional dependencies from
`requirements.txt`. Historical environment records are in `docs/`.
The local network-mesh helper's attribution is in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md);
no upstream patch or `networks_fenicsx` installation is required.

Set these before submission to override paths:

* `STUDY_VENV`: defaults to `/projects/br1/pzz1/segmentation-sensitivity/software/venv`.
* `STUDY_WORK`: NOTS scratch base, default `/scratch/pzz1/segmentation-sensitivity`. Each experiment uses a named subdirectory.
* `STUDY_DATA_ZIP`: defaults to `/projects/br1/pzz1/segmentation-sensitivity/oilblooddata.zip`.

Without the NOTS environment or `STUDY_WORK`, experiment scripts put raw work in
`experiments/<name>/results/raw/`. Published outputs always stay with their
experiment. The library has no experiment or repository-path dependencies.

## Existing results

To refresh published reports from the earlier completed results without rerunning
simulations:

```bash
source environment.sh
python -m experiments.collect --from-results /home/pzz1/segmentation-sensitivity-results
```

No workflow creates a tarball. To copy the repository to a local computer:

```bash
rsync -av --exclude='.work/' --exclude='__pycache__/' --exclude='.git/' \
  pzz1@nots.rice.edu:/home/pzz1/segmentation-sensitivity/ ./segmentation-sensitivity/
```

All potentials and fluxes are dimensionless; VTU coordinates are millimetres.
The original ZIP and historical full-resolution raw solutions remain on NOTS.
The recorded numerical regression predates the submission-interface refactor;
changing the submitter does not rerun simulations.

Reference: [Laurino and Zunino (2019)](https://doi.org/10.1051/m2an/2019042),
equations (11)–(12). Parameters, validation, and limitations are in the reports.
