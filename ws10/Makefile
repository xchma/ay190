TEXS   = ws10.tex
PLOT   = fig1.pdf fig2.pdf
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
