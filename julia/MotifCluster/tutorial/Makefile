all: tutorial

.PHONY: tutorial

tutorial:
	lualatex -pdf motifcluster_tutorial.tex
	bibtex motifcluster_tutorial
	pythontex motifcluster_tutorial.tex
	lualatex -pdf motifcluster_tutorial.tex
	lualatex -pdf motifcluster_tutorial.tex
