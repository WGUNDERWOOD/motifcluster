# motif-based-clustering

[![Build Status](https://travis-ci.com/WGUNDERWOOD/motif-based-clustering.svg?branch=master)](https://travis-ci.com/github/WGUNDERWOOD/motif-based-clustering)
[![Python Docs Status](https://img.shields.io/readthedocs/motifcluster?label=python%20docs)](https://motifcluster.readthedocs.io/en/latest/)
[![Python Coverage Status](https://img.shields.io/coveralls/github/WGUNDERWOOD/motif-based-clustering?label=python%20coverage)](https://coveralls.io/github/WGUNDERWOOD/motif-based-clustering)
[![R Coverage Status](https://img.shields.io/codecov/c/github/wgunderwood/motif-based-clustering?label=R%20coverage)](https://codecov.io/gh/WGUNDERWOOD/motif-based-clustering)

Source code for my paper
*Motif-Based Spectral Clustering of Weighted Directed Networks*,
with Andrew Elliott and Mihai Cucuringu.
A preprint of this paper is avaliable at
[arXiv:2004.01293](https://arxiv.org/abs/2004.01293).

## Introduction

This repository provides implementations of motif-based spectral clustering
of weighted directed networks in R and in Python.
These provide the capability for:

- Building motif adjacency matrices
- Sampling random weighted directed networks
- Spectral embedding with motif adjacency matrices
- Motif-based spectral clustering

The methods are all designed to run quickly on large sparse networks,
and are easy to install and use.

### Branches

The master branch contains stable versions.
The develop branch may be unstable,
and is for development purposes only.

### License

This repository,
along with its included R package and Python package,
are all licensed under
[GPL-3](http://gplv3.fsf.org/).





## R package

The R package "motifcluster" is in the [R](./R/) directory.

### Installation

The R package can be installed from the GitHub master branch with:

```
install_github("wgunderwood/motif-based-clustering/R")
```

### Dependencies

The R package has the following dependencies, available on CRAN:

- igraph
- LICORS
- Matrix
- RSpectra

### Documentation

The package's manual and an instructional vignette are in the
[R/doc](./R/doc) directory.
R documentation files are provided for each function
available in the package.




## Python code

The Python package "motifcluster" is in the
[python](./python/) directory.

### Installation

The Python package can be installed from PyPI with

```
pip install motifcluster
```

Or, to install to the user directory, with

```
pip install --user motifcluster
```

### Dependencies

The Python package has the following dependencies,
available on PyPI:

- Networkx
- Numpy
- Scipy
- Scikit-learn

### Documentation

The package's manual is in the
[python/doc](./python/doc/) directory.
Full documentation is also available at
[motifcluster.readthedocs.io](https://motifcluster.readthedocs.io/).
