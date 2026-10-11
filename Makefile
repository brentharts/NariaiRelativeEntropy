# NariaiRelativeEntropy -- build system
#
#   make              -> the paper (PAPER below), built with pdflatex twice
#   make paper        -> same as above
#   make supplementary -> supplementary.pdf (figures + TikZ document)
#   make figures      -> regenerate matplotlib figures into figures/
#   make verify       -> run the Appendix B numerical checks
#   make selftest     -> pass/fail assertions on every claim
#   make convergence  -> grid and truncation study
#   make lean         -> regenerate and kernel-check NariaiFacts.lean
#   make exact        -> exact-spectrum packet + Fock relative entropy
#   make check        -> the referee-report checks and the open-item checks
#   make clean        -> remove build artifacts (keeps PDFs)
#   make distclean    -> remove build artifacts and PDFs

PY      ?= python3
LATEX   ?= pdflatex -interaction=nonstopmode -halt-on-error

FIGS    := figures/fig_flux_windows.pdf \
           figures/fig_dichotomy.pdf \
           figures/fig_one_mode.pdf

# Current manuscript: change this one line when a new version is written.
PAPER   ?= nariai_relative_entropy_v14

OPEN    := $(shell command -v xdg-open 2>/dev/null || command -v open 2>/dev/null)

default: paper

paper:
	$(LATEX) $(PAPER).tex
	$(LATEX) $(PAPER).tex

view: paper
	$(OPEN) $(PAPER).pdf

supplementary:
	$(LATEX) supplementary.tex
	$(LATEX) supplementary.tex


# --- figures -----------------------------------------------------------
$(FIGS) &: generate_figures.py
	$(PY) generate_figures.py

figures: $(FIGS)

# --- numerical verification (the checks the paper refers to) --------------
verify: appendix_b_verification.py
	$(PY) appendix_b_verification.py

selftest: appendix_b_verification.py
	$(PY) appendix_b_verification.py --selftest

convergence: appendix_b_verification.py
	$(PY) appendix_b_verification.py --convergence

# --- machine-checked fragments -----------------------------------------
lean: nariailean.py nariai_exact.py
	$(PY) nariailean.py --check

exact: nariai_exact.py
	$(PY) nariai_exact.py

check: check_referee_points.py check_open_items.py check_normalization.py
	$(PY) check_referee_points.py
	$(PY) check_open_items.py
	$(PY) check_normalization.py





# --- housekeeping --------------------------------------------------------
clean:
	rm -f *.aux *.log *.out *.toc *.synctex.gz

distclean: clean
	rm -f supplementary.pdf $(PAPER).pdf
	rm -rf figures/ __pycache__/

.PHONY: default paper view supplementary figures verify selftest convergence lean exact check clean distclean
