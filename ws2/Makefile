TEXS   = ws2.tex
TABLE  = tb1.tex
PLOTS  = fig1.pdf fig2.pdf
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
