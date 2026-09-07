# Build the report. `make` produces report.pdf from the committed figures.
#
#   make          - build report.pdf
#   make figures  - regenerate the figures from figures/data/*.csv
#   make data     - regenerate that data by re-running the solvers (hours)
#   make clean    - remove LaTeX intermediates
#
# The measurements are committed under figures/data so that the report builds
# without a recompute; `make data` is how they were produced.

PYTHON ?= python3
ROOT := ../..
FIGURES := $(patsubst %,figures/fig_%.pdf,sweep history profiles currents recolouring)

.PHONY: all figures data clean distclean

all: report.pdf

report.pdf: report.tex $(FIGURES)
	latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

figures $(FIGURES): figures/make_figures.py $(wildcard figures/data/*.csv)
	$(PYTHON) figures/make_figures.py

data: report_data
	./report_data sweep    figures/data
	./report_data history  figures/data
	./report_data profiles figures/data
	./report_data currents figures/data

report_data: tools/report_data.cpp
	g++ -O2 -std=c++17 -ffp-contract=off -fopenmp -I $(ROOT)/src -o $@ $< \
	    $(ROOT)/src/lbm/*.cpp $(ROOT)/src/omp/*.cpp

clean:
	latexmk -C report.tex 2>/dev/null || true
	rm -f report.aux report.log report.out report.toc report.fls report.fdb_latexmk

distclean: clean
	rm -f report.pdf report_data $(FIGURES)
