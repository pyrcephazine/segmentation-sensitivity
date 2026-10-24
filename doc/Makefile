# Build the methods note. Requires a TeX Live install with pdflatex,
# and (only if regenerating the figure) python3 with numpy + matplotlib.

TEX = vessel_break_equivalence
FIG = phi_figure.pdf

.PHONY: all figure clean

# Default: build the PDF from the existing figure (phi_figure.pdf is shipped).
all: $(TEX).pdf

$(TEX).pdf: $(TEX).tex $(FIG)
	pdflatex -interaction=nonstopmode -halt-on-error $(TEX).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(TEX).tex

# Regenerate the two-panel phi figure from source (optional).
figure: make_fig.py
	python3 make_fig.py

# Recompute the Table 1 / worked-example numbers (prints to stdout).
numbers: phi_demo.py
	python3 phi_demo.py

clean:
	rm -f $(TEX).aux $(TEX).log $(TEX).out $(TEX).pdf
