TEXS   = bibtex_example.tex
BIB    = $(TEXS:.tex=.bib)
SOURCE = $(TEXS) $(BIB)
OUTP   = $(TEXS:.tex=.pdf)
CRAP   = $(TEXS:.tex=.aux) $(TEXS:.tex=.log) $(TEXS:.tex=.bbl) $(TEXS:.tex=.blg)
COMP   = pdflatex
COMPB  = bibtex

$(OUTP): $(SOURCE)
	$(COMP) $(TEXS)
	$(COMPB) $(BIB:.bib=)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)

clean:
	rm $(CRAP)
