.libPaths("~/R_packages")
Sys.setenv(RSTUDIO_PANDOC="/usr/bin/pandoc")
library(devtools)
library(roxygen2)
library(lintr)
library(covr)

#check()
#build_manual(path = "./doc")
#build_vignettes()
report()
#lint_package(linters = with_defaults(
               #object_name_linter = NULL,
               #object_usage_linter = NULL,
               #cyclocomp_linter = NULL))
