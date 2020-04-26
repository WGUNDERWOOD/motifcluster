#!/usr/bin/env bash

# check running in right directory
if [ ${PWD##*/} != "motif-based-clustering" ]
then
  echo "wrong directory!"
  exit
fi

# base directory
basedir=$PWD

# python
echo "Python"
cd $basedir/python/

# python testing
pytest --profile-svg

# python linting
pylint_ignore=_build_G,_build_J,_build_Gs,_build_Js,
pylint_ignore+=_build_Gd,_build_Jd,_build_J0,_build_Jn,
pylint_ignore+=_build_Id,_build_Je,_build_Gp,
pylint_ignore+=_a_one_b,_a_b_one

pylint \
  --indent-string="  " \
  --function-naming-style=any \
  --variable-naming-style=any \
  --exclude-protected=$pylint_ignore \
  motifcluster

# python distribution packages
python setup.py sdist bdist_wheel

# python rtfd
cd $basedir/python/doc/
make html latex

# R
echo "R"
cd $basedir/R/
Rscript build_R.R
