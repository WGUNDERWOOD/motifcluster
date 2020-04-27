#!/usr/bin/env bash

# python testing
pytest -n 8 --profile-svg

# python linting
pylint --rcfile=.pylintrc motifcluster

# python distribution packages
python setup.py sdist bdist_wheel

# python rtfd
cd doc/
make html latex latexpdf
cp _build/latex/motifcluster.pdf .
