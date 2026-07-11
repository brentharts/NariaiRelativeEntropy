# NariaiRelativeEntropy -- build system
#
#   make              -> supplementary.pdf (figures + TikZ document)
#   make figures      -> regenerate matplotlib figures into figures/
#   make verify       -> run the Appendix B numerical checks
#   make docs     -> supplementary
#   make clean        -> remove build artifacts (keeps PDFs)
#   make distclean    -> remove build artifacts and PDFs

PY      ?= python3
LATEX   ?= pdflatex -interaction=nonstopmode -halt-on-error

FIGS    := figures/fig_flux_windows.pdf \
           figures/fig_dichotomy.pdf \
           figures/fig_one_mode.pdf

default:
	pdflatex supplementary.tex
	pdflatex supplementary.tex
	open supplementary.pdf


# --- figures -----------------------------------------------------------
$(FIGS) &: generate_figures.py
	$(PY) generate_figures.py

figures: $(FIGS)

# --- numerical verification (Appendix B of the main paper) --------------
verify: appendix_b_verification.py
	$(PY) appendix_b_verification.py





# --- housekeeping --------------------------------------------------------
clean:
	rm -f *.aux *.log *.out *.toc *.synctex.gz

distclean: clean
	rm -f supplementary.pdf
	rm -rf figures/ __pycache__/
