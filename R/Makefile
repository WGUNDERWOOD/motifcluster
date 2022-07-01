# Makefile: R

all: doc check man vignette cov lint

.PHONY: setup doc check man vignette cov lint install release clean

setup:
	@echo ".libPaths('~/R_packages');" > .setup.R
	@echo "Sys.setenv(RSTUDIO_PANDOC='/usr/bin/pandoc');" >> .setup.R
	@echo "library(devtools);" >> .setup.R
	@echo "library(roxygen2);" >> .setup.R
	@echo "library(lintr);" >> .setup.R
	@echo "library(covr);" >> .setup.R
	@echo "library(filesstrings);" >> .setup.R

doc: setup
	@cat .setup.R > .temp.R
	@echo "document()" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

check: setup
	@cat .setup.R > .temp.R
	@echo "check()" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

man: setup
	@cat .setup.R > .temp.R
	@echo "build_manual(path = './doc')" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

vignette: setup
	@cat .setup.R > .temp.R
	@echo "build_vignettes()" >> .temp.R
	@Rscript .temp.R
	@rm doc/motifcluster_vignette.R doc/motifcluster_vignette.Rmd
	@mv doc/motifcluster_vignette.pdf vignettes/
	@rm .gitignore
	@rm .setup.R .temp.R

cov: setup
	@cat .setup.R > .temp.R
	@echo "cov <- package_coverage(); zero_coverage(cov)" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

lint: setup
	@cat .setup.R > .temp.R
	@echo "lint_package(linters = with_defaults( \
		object_name_linter = NULL, \
		object_usage_linter = NULL, \
		cyclocomp_linter = NULL))" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

install: setup
	@cat .setup.R > .temp.R
	@echo "install()" >> .temp.R
	@Rscript .temp.R
	@rm .setup.R .temp.R

release: setup
	@cat .setup.R > .temp.R
	@echo "release()" >> .temp.R
	@R -f .temp.R
	@rm .setup.R .temp.R

clean:
	@rm -f .setup.R .temp.R
	@rm -f .gitignore
	@rm -rf doc/
	@rm -rf man/
	@rm -rf Meta/
