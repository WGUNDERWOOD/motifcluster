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

# python tutorial
echo "Building python tutorial"
cd ../tutorial
bash build_tutorial.sh

# python distribution packages
echo "Building python distribution packages"
cd ..
python setup.py sdist bdist_wheel
python -m twine check dist/*

# pip uninstall
#pip uninstall motifcluster

# python install local
#pip install --user dist/motifcluster-0.1.2.tar.gz


# python upload to PyPI
#python -m twine upload --repository testpypi dist/*
#python -m twine upload dist/*
### username: __token__
### password: <PyPI API token>
