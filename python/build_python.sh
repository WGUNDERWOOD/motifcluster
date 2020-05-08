#!/usr/bin/env bash

echo "Running python tests"
pytest --cov=motifcluster --profile-svg tests/

# python linting
echo "Running python linter"
pylint -j 8 --rcfile=.pylintrc motifcluster

# python rtfd
echo "Building python docs"
cd doc/
make html

# python distribution packages
cd ..
echo "Building python distribution packages"
python setup.py sdist bdist_wheel
python -m twine check dist/*

# python upload to PyPI
#python -m twine upload --repository testpypi dist/*
#python -m twine upload dist/*
### username: __token__
### password: <PyPI API token>
