[![build status](https://travis-ci.com/WGUNDERWOOD/motifcluster.svg?branch=master)](https://travis-ci.com/github/WGUNDERWOOD/motifcluster)
[![license: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/249007000.svg)](https://zenodo.org/badge/latestdoi/249007000)
[![CRAN version](https://img.shields.io/cran/v/motifcluster?color=7733BB&label=CRAN)](https://cran.r-project.org/web/packages/motifcluster/index.html)
[![PyPI version](https://img.shields.io/pypi/v/motifcluster?color=7733BB&label=PyPI)](https://pypi.org/project/motifcluster/)
[![R coverage status](https://img.shields.io/codecov/c/github/wgunderwood/motifcluster?label=R%20cov)](https://codecov.io/gh/WGUNDERWOOD/motifcluster)
[![python coverage status](https://img.shields.io/coveralls/github/WGUNDERWOOD/motifcluster?label=python%20cov)](https://coveralls.io/github/WGUNDERWOOD/motifcluster)
[![python docs status](https://img.shields.io/readthedocs/motifcluster?label=python%20docs)](https://motifcluster.readthedocs.io/en/latest/)


# motifcluster <img src="https://github.com/WGUNDERWOOD/motifcluster/raw/develop/sticker/hex_sticker_small.png" alt="motifcluster sticker" align="right" />



Motif-based spectral clustering of weighted directed networks

## Introduction

This repository provides implementations of motif-based spectral clustering
of weighted directed networks in R and in Python.
This code is based on methods detailed in
[Underwood, Elliott and Cucuringu, 2020],
which is available at
[arXiv:2004.01293](https://arxiv.org/abs/2004.01293).
These packages provide the capability for:

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

### Authors

  - [William George Underwood](https://wgunderwood.github.io/),
    Princeton University
    (Python, R, maintainer)
  - [Andrew Elliott](https://www.turing.ac.uk/people/researchers/andrew-elliott),
    The Alan Turing Institute
    (Python)

### License

This repository,
along with its included R package and Python package,
are all licensed under
[GPLv3](http://gplv3.fsf.org/).





## R package

The **motifcluster** R package is in the
[R](https://github.com/WGUNDERWOOD/motifcluster/tree/master/R)
directory.

### Installation

The R package can be installed from CRAN with:

```
install.packages("motifcluster")
```

### Dependencies

The R package has the following dependencies, available on CRAN:

- igraph
- LICORS
- Matrix
- RSpectra

### Documentation

The package's manual is in the
[R/doc](https://github.com/WGUNDERWOOD/motifcluster/tree/master/R/doc)
directory.
R documentation files are provided for each function
available in the package.
An instructional vignette is in the
[R/vignettes](https://github.com/WGUNDERWOOD/motifcluster/tree/master/R/vignettes)
directory.




## Python package

The **motifcluster** Python package is in the
[python](https://github.com/WGUNDERWOOD/motifcluster/tree/master/python)
directory.

### Installation

The Python package can be installed from PyPI with:

```
pip install motifcluster
```

### Dependencies

The Python package has the following dependencies,
available on PyPI:

- Networkx
- Numpy
- Scipy
- Scikit-learn

### Documentation

Full documentation is available at
[motifcluster.readthedocs.io](https://motifcluster.readthedocs.io/).

## Sticker

A high-resolution hexagonal sticker is available in the
[sticker](https://github.com/WGUNDERWOOD/motifcluster/tree/master/sticker)
directory.
