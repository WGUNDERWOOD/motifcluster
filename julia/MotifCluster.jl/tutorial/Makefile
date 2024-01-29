all: tutorial

.PHONY: tutorial clean

tutorial:
	latexmk -pdf -quiet motifcluster_tutorial.tex
	pythontex motifcluster_tutorial.tex
	latexmk -pdf -quiet motifcluster_tutorial.tex

clean:
	rm -rfv \
		auto/ \
		*.aux \
		*-blx.bib \
		*.bbl \
		*.bcf \
		*.blg \
		*.dvi \
		*.fdb_latexmk \
		*.fls \
		*.glg \
		*.glo \
		*.gls \
		*.ist \
		*.lof \
		*.log \
		*.lot \
		*.nav \
		*.out \
		pythontex-files-*/ \
		*.pytxcode \
		*.run.xml \
		*.sha256 \
		*.snm \
		*.synctex.gz \
		*.toc \
		*.xdv \
