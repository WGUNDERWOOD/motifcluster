#!/usr/bin/env bash

# check running in right directory
if [ ${PWD##*/} != "motifcluster" ]
then
  echo "Running build.sh from the wrong directory!"
  exit
fi

# base directory
basedir=$PWD

# python
echo "Python build"
cd $basedir/python/
bash build_python.sh

# R
echo "R build"
cd $basedir/R/
Rscript build_R.R

# Compress pdf files
echo "Compress pdfs"
cd $basedir
bash compress_pdfs.sh

# Performance
echo "Testing performance"
cd performance/
bash performance.sh

# Finished
echo "Build finished!"
