# motif-based-clustering

[![Build Status](https://img.shields.io/travis/wgunderwood/motif-based-clustering/master.svg?label=master)](https://travis-ci.com/github/WGUNDERWOOD/motif-based-clustering)
[![Build Status](https://img.shields.io/travis/wgunderwood/motif-based-clustering/develop.svg?label=develop)](https://travis-ci.com/github/WGUNDERWOOD/motif-based-clustering)

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



## R package

The R package `motifcluster` is in the [R](./R/) directory.

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

The Python source code is in the [python](./python/) directory
and requires the following Python packages:

- Networkx
- Numpy
- Scipy
