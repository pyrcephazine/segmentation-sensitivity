# pt002 inflow

One intact supplied inflow network coupled to its supplied tissue mesh, plus a
vessel-resolution check. `run.py` imports the ZIP members, normalizes coordinates,
selects the positive inlet, and runs the shared solver. `report.py` exports VTK
solutions, checks conservation and refinement, and writes the report.

From the repository root after sourcing `environment.sh`:

```bash
python -m experiments.pt002_inflow.run --prepare
python3 submit.py supplied
```

Run the synthetic validation first: the supplied report includes a comparison
against its saved line-case baseline. See the root README for the full workflow.
Published results are in `results/`, including `REPORT.html`, `REPORT.pdf`,
`metrics.csv`, and `coupled_solution.vtu`. With `STUDY_WORK` set, raw work is under
`$STUDY_WORK/pt002_inflow/`; otherwise it is in `results/raw/`.
