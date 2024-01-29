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

man: .setup.R
	@echo -e "\e[0;35m\033[1mMaking R man pages...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "document()" >> .temp.R
	@Rscript .temp.R

check: .setup.R
	@echo -e "\e[0;35m\033[1mChecking R package...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "check()" >> .temp.R
	@Rscript .temp.R

test: .setup.R
	@echo -e "\e[0;35m\033[1mTesting R package...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "test()" >> .temp.R
	@Rscript .temp.R

doc: .setup.R man
	@echo -e "\e[0;35m\033[1mMaking R documentation...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@mkdir -p doc/
	@echo "build_manual(path = './doc')" >> .temp.R
	@Rscript .temp.R
	@echo

vignette: .setup.R
	@echo -e "\e[0;35m\033[1mMaking R vignette...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "build_vignettes()" >> .temp.R
	@Rscript .temp.R
	@rm doc/motifcluster_vignette.R doc/motifcluster_vignette.Rmd
	@mv doc/motifcluster_vignette.pdf vignettes/
	@rm .gitignore

cov: .setup.R
	@echo -e "\e[0;35m\033[1mChecking R coverage...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "cov <- package_coverage(); zero_coverage(cov)" >> .temp.R
	@Rscript .temp.R

lint: .setup.R
	@echo -e "\e[0;35m\033[1mLinting R package...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "lint_package(linters = linters_with_defaults( \
		object_name_linter = NULL, \
		object_usage_linter = NULL, \
		indentation_linter = NULL, \
		cyclocomp_linter = NULL))" >> .temp.R
	@Rscript .temp.R

install: .setup.R
	@echo -e "\e[0;35m\033[1mInstalling R package...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "install()" >> .temp.R
	@Rscript .temp.R

release: .setup.R
	@echo -e "\e[0;35m\033[1mReleasing R package...\e[0;30m\033[0m"
	@cat .setup.R > .temp.R
	@echo "release()" >> .temp.R
	@R -f .temp.R

clean:
	@echo -e "\e[0;35m\033[1mCleaning up R package...\e[0;30m\033[0m"
	@rm -f .setup.R .temp.R
	@rm -f .gitignore
	@rm -rf man/
	@rm -rf Meta/
