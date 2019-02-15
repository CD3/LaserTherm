BASE=$(shell basename $$PWD )

$(BASE).pdf: main.dvi
	dvipdf main.dvi main.pdf
	gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=$@ main.pdf

main.dvi: main.tex
	latex main.tex
	latex main.tex

