#!/usr/bin/env bash

# check running in right directory
if [ ${PWD##*/} != "motifcluster" ]
then
  echo "wrong directory!"
  exit
fi

# base directory
basedir=$PWD

# python
echo "Python"
cd $basedir/python/
bash build_python.sh | tee ../build.log

# R
echo "R"
cd $basedir/R/
Rscript build_R.R | tee -a ../build.log

# Compress pdf files
echo "Compress pdfs"
cd $basedir
bash compress_pdfs.sh | tee -a build.log

# Finished
printf "\n\n\n"
echo "Build finished! Looking for errors..."
printf "\n\n\n"
grep -iC 10 "error\|warn\|note" build.log
