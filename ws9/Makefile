TEXS   = ws9.tex
TABLE  = tb.tex
PLOT   = fig.png
SOURCE = $(TEXS) $(TABLE) $(PLOT)
OUTP   = $(TEXS:.tex=.pdf)
CRAP   = $(TEXS:.tex=.aux) $(TEXS:.tex=.log)
COMP   = pdflatex

$(OUTP): $(SOURCE)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)

clean:
	rm $(CRAP)
