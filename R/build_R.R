.libPaths("~/R_packages")
Sys.setenv(RSTUDIO_PANDOC="/usr/bin/pandoc")
library(devtools)
library(roxygen2)
library(lintr)
library(covr)
library(filesstrings)

print("Building R docs")
document()

print("Checking R package")
check()

print("Building R manual")
build_manual(path = "./doc")

print("Building R vignettes")
build_vignettes()
file.remove("doc/motifcluster_vignette.R")
file.remove("doc/motifcluster_vignette.Rmd")
file.move("doc/motifcluster_vignette.pdf", "vignettes")

print("Checking R coverage")
cov <- package_coverage()
zero_coverage(cov)

print("Running R linter")
lint_package(linters = with_defaults(
               object_name_linter = NULL,
               object_usage_linter = NULL,
               cyclocomp_linter = NULL))

print("Removing R gitignore")
file.remove(".gitignore")

#print("Installing package")
#install()

#print("Releasing")
#release()
