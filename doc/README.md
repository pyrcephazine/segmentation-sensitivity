# Vessel-break flow-equivalence methods note

Schur-complement port reduction of a 3D-bridged vessel break and the
resistance-fraction (phi) tolerance, sitting next to the Amare et al. 1D
Hagen-Poiseuille formulation.

## Build

    make            # -> vessel_break_equivalence.pdf  (runs pdflatex twice)

Requires TeX Live with `pdflatex` and the TikZ libraries used in the preamble
(arrows.meta, positioning, calc, backgrounds, fit, shapes.geometric,
shapes.misc) plus amsmath/amsthm, bm, booktabs, graphicx, hyperref — all
standard in a full TeX Live.

## Files

- `vessel_break_equivalence.tex` — the note. Figure 1 (break geometry, tet-mesh
  ports/interior) and Figure 2 (Table-1 break locations on the path) are inline
  TikZ; Figure 3 is the external `phi_figure.pdf`.
- `phi_figure.pdf` — two-panel figure (flow-ratio vs. bridge error; phi by
  generation). Shipped prebuilt so `make` needs no Python.
- `make_fig.py` — regenerates `phi_figure.pdf`  (`make figure`; needs numpy,
  matplotlib).
- `phi_demo.py` — recomputes the worked-example numbers in Table 1
  (`make numbers`; needs numpy).

## Regenerating everything from scratch

    make figure     # rebuild phi_figure.pdf from make_fig.py
    make            # rebuild the PDF
