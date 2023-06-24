all: tutorial

.PHONY: tutorial

tutorial:
	lualatex -pdf motifcluster_tutorial.tex
	bibtex motifcluster_tutorial
	pythontex motifcluster_tutorial.tex
	lualatex -pdf motifcluster_tutorial.tex
	lualatex -pdf motifcluster_tutorial.tex

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
