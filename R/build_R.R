.libPaths("~/R_packages")
Sys.setenv(RSTUDIO_PANDOC="/usr/bin/pandoc")
library(devtools)
library(roxygen2)
library(lintr)
library(covr)

<<<<<<< Updated upstream
#check()
#build_manual(path = "./doc")
#build_vignettes()
cov <- package_coverage()
zero_coverage(cov)
#lint_package(linters = with_defaults(
               #object_name_linter = NULL,
               #object_usage_linter = NULL,
               #cyclocomp_linter = NULL))
=======
check()
build_manual(path = "./doc")
build_vignettes()
cov <- package_coverage()
zero_coverage(cov)
lint_package(linters = with_defaults(
               object_name_linter = NULL,
               object_usage_linter = NULL,
               cyclocomp_linter = NULL))
>>>>>>> Stashed changes
