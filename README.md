[![CI](https://github.com/WGUNDERWOOD/motifcluster/actions/workflows/CI.yml/badge.svg)](https://github.com/WGUNDERWOOD/motifcluster/actions/workflows/CI.yml)
[![license: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](https://img.shields.io/cran/v/motifcluster?color=7733BB&label=CRAN)](https://cran.r-project.org/web/packages/motifcluster/index.html)
[![PyPI version](https://img.shields.io/pypi/v/motifcluster?color=7733BB&label=PyPI)](https://pypi.org/project/motifcluster/)
[![codecov](https://codecov.io/gh/WGUNDERWOOD/motifcluster/branch/main/graph/badge.svg?token=DbGSOsocw6)](https://codecov.io/gh/WGUNDERWOOD/motifcluster)
[![python docs status](https://img.shields.io/readthedocs/motifcluster?label=python%20docs)](https://motifcluster.readthedocs.io/en/latest/)
[![Julia build Status](https://github.com/WGUNDERWOOD/DyadicKDE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WGUNDERWOOD/DyadicKDE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Julia docs status](https://img.shields.io/badge/docs-stable-blue.svg)](https://WGUNDERWOOD.github.io/DyadicKDE.jl/stable)
[![Julia Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

TODO split badges by language






# motifcluster <img src="https://github.com/WGUNDERWOOD/motifcluster/raw/main/sticker/hex_sticker_small.png" alt="motifcluster sticker" align="right" width=120 />



Motif-based spectral clustering of weighted directed networks

## Introduction

This repository provides implementations of motif-based spectral clustering
of weighted directed networks in R, Python and Julia.
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

The main branch contains stable versions.
The develop branch may be unstable,
and is for development purposes only.

### Authors

  - [William George Underwood](https://wgunderwood.github.io/),
    Princeton University
    (Python, R, Julia, maintainer)
  - [Andrew Elliott](https://www.turing.ac.uk/people/researchers/andrew-elliott),
    The Alan Turing Institute
    (Python)

### License

This repository
and its included R, Python and Julia packages
are all licensed under
[GPLv3](http://gplv3.fsf.org/).





## R package

The **motifcluster** R package is in the
[R](https://github.com/WGUNDERWOOD/motifcluster/tree/main/R)
directory.

### Installation

The R package can be installed from CRAN with:

```
install.packages("motifcluster")
```

### Dependencies

The R package has the following dependencies, available on CRAN:

- igraph
- Matrix
- RSpectra

### Documentation

The package's manual is in the
[R/doc](https://github.com/WGUNDERWOOD/motifcluster/tree/main/R/doc)
directory.
R documentation files are provided for each function
available in the package.
An instructional vignette is in the
[R/vignettes](https://github.com/WGUNDERWOOD/motifcluster/tree/main/R/vignettes)
directory.




## Python package

The **motifcluster** Python package is in the
[python](https://github.com/WGUNDERWOOD/motifcluster/tree/main/python)
directory.

### Installation

The Python package can be installed from PyPI with:

```
pip install motifcluster
```

Alternatively it can be installed with conda using:

```
conda install -c conda-forge motifcluster
```


### Dependencies

The Python package has the following dependencies:

- Networkx
- Numpy
- Scipy
- Scikit-learn

### Documentation

Full documentation is available at
[motifcluster.readthedocs.io](https://motifcluster.readthedocs.io/).




## Julia package

The **MotifCluster** Julia package is in the
[julia](https://github.com/WGUNDERWOOD/motifcluster/tree/main/julia)
directory.

### Installation

From the Julia General registry:

```
] add MotifCluster
```

### Dependencies
- Aqua
- Clustering
- Distributions
- Graphs

### Documentation
Documentation for the **MotifCluster** package is available on 
[the web](https://wgunderwood.github.io/MotifCluster.jl/stable/).







## Performance

The
[performance](https://github.com/WGUNDERWOOD/motifcluster/tree/main/performance)
directory contains scripts and plots relating to timing
the construction of motif adjacency matrices,
in both R and Python.

## Sticker

A high-resolution hexagonal sticker is available in the
[sticker](https://github.com/WGUNDERWOOD/motifcluster/tree/main/sticker)
directory.
