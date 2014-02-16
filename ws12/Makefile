TEXS   = ws12.tex
PLOT   = fig*
SOURCE = $(TEXS) $(PLOT)
OUTP   = $(TEXS:.tex=.pdf)
CRAP   = $(TEXS:.tex=.aux) $(TEXS:.tex=.log)
COMP   = pdflatex

$(OUTP): $(SOURCE)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)

clean:
	rm $(CRAP)
