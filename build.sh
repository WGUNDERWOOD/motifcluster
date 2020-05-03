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
bash build_python.sh

# R
echo "R"
cd $basedir/R/
Rscript build_R.R
rm .gitignore

# Compress pdf files
echo "Compress pdfs"
cd $basedir
bash compress_pdfs.sh
