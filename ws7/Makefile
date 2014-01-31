TEXS   = ws7.tex
PLOTS  = fig1.pdf fig2.pdf
TABLE  = tb1.tex tb2.tex tb3.tex
SOURCE = $(TEXS) $(TABLE) $(PLOTS)
OUTP   = $(TEXS:.tex=.pdf)
CRAP   = $(TEXS:.tex=.aux) $(TEXS:.tex=.log)
COMP   = pdflatex

$(OUTP): $(SOURCE)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)
	$(COMP) $(TEXS)

clean:
	rm $(CRAP)
