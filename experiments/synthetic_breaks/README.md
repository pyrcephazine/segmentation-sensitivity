# Synthetic breaks

Six line cases, five Y cases, and 31 tree cases (intact plus five seeds at six
nested break counts). `cases.py` defines this experiment's protocol; `run.py`
uses the reusable `segmentation_sensitivity` library.

From the repository root after sourcing `environment.sh`, run `python3 submit.py synthetic`.
The submitter validation jobs reproduce the analytic, MPI, and refinement
checks. See the repository README for the complete sequence.

Published outputs are in `results/REPORT.html`, `results/REPORT.pdf`, and adjacent
figures and tables. `config.py` sets the results and raw-work locations. With
`STUDY_WORK` set, raw work is under `$STUDY_WORK/synthetic_breaks/`; otherwise it
is in `results/raw/`.
