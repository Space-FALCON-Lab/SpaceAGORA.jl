# Makefile for AIAA paper on RL-based constellation seeding

LATEX = pdflatex
BIBTEX = bibtex
MAIN = main

.PHONY: all clean view

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex introduction.tex background.tex methodology.tex embeddings.tex results.tex conclusion.tex acknowledgments.tex references.bib
	$(LATEX) $(MAIN)
	$(BIBTEX) $(MAIN)
	$(LATEX) $(MAIN)
	$(LATEX) $(MAIN)

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.pdf

view: $(MAIN).pdf
	xdg-open $(MAIN).pdf 2>/dev/null || open $(MAIN).pdf 2>/dev/null || start $(MAIN).pdf 2>/dev/null
