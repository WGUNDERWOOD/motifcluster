# Makefile: R

all: doc test check man vignette cov lint

.PHONY: doc test check man vignette cov lint install release clean

.setup.R:
	@echo ".libPaths('~/R_packages');" > .setup.R
	@echo "Sys.setenv(RSTUDIO_PANDOC='/usr/bin/pandoc');" >> .setup.R
	@echo "library(devtools);" >> .setup.R
	@echo "library(roxygen2);" >> .setup.R
	@echo "library(lintr);" >> .setup.R
	@echo "library(covr);" >> .setup.R
	@echo "library(filesstrings);" >> .setup.R

doc: .setup.R
	@cat .setup.R > .temp.R
	@echo "document()" >> .temp.R
	@Rscript .temp.R

check: .setup.R
	@cat .setup.R > .temp.R
	@echo "check()" >> .temp.R
	@Rscript .temp.R

test: .setup.R
	@cat .setup.R > .temp.R
	@echo "test()" >> .temp.R
	@Rscript .temp.R

man: .setup.R doc
	@cat .setup.R > .temp.R
	@mkdir -p doc/
	@echo "build_manual(path = './doc')" >> .temp.R
	@Rscript .temp.R
	@echo
	@for f in doc/motifcluster_0.2.0.pdf ; do optpdf $$f; optpdf $$f ; done

vignette: .setup.R
	@cat .setup.R > .temp.R
	@echo "build_vignettes()" >> .temp.R
	@Rscript .temp.R
	@rm doc/motifcluster_vignette.R doc/motifcluster_vignette.Rmd
	@optpdf doc/motifcluster_vignette.pdf
	@optpdf doc/motifcluster_vignette.pdf
	@mv doc/motifcluster_vignette.pdf vignettes/
	@rm .gitignore

cov: .setup.R
	@cat .setup.R > .temp.R
	@echo "cov <- package_coverage(); zero_coverage(cov)" >> .temp.R
	@Rscript .temp.R

lint: .setup.R
	@cat .setup.R > .temp.R
	@echo "lint_package(linters = with_defaults( \
		object_name_linter = NULL, \
		object_usage_linter = NULL, \
		cyclocomp_linter = NULL))" >> .temp.R
	@Rscript .temp.R

install: .setup.R
	@cat .setup.R > .temp.R
	@echo "install()" >> .temp.R
	@Rscript .temp.R

release: .setup.R
	@cat .setup.R > .temp.R
	@echo "release()" >> .temp.R
	@R -f .temp.R

clean:
	@rm -f .setup.R .temp.R
	@rm -f .gitignore
	@rm -rf doc/
	@rm -rf man/
	@rm -rf Meta/
