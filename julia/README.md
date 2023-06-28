# motifcluster <img src="https://github.com/WGUNDERWOOD/motifcluster/raw/main/sticker/hex_sticker_small.png" alt="motifcluster sticker" align="right" width=120 />

Motif-based spectral clustering of weighted directed networks in Julia.

## Introduction

The **MotifCluster** package provides
implementations of motif-based spectral clustering
of weighted directed networks in Julia.
These provide the capability for:

- Building motif adjacency matrices
- Sampling random weighted directed networks
- Spectral embedding with motif adjacency matrices
- Motif-based spectral clustering

The methods are all designed to run quickly on large sparse networks,
and are easy to install and use.
These methods are based on those described in
[Underwood, Elliott and Cucuringu, 2020],
which is available at
[arXiv:2004.01293](https://arxiv.org/abs/2004.01293).

## Installation

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


## Tutorial

An instructional tutorial for the **MotifCluster** package
is available in the
[tutorial](https://github.com/WGUNDERWOOD/motifcluster/tree/main/julia/MotifCluster/tutorial)
directory.

## Author

- [William George Underwood](https://wgunderwood.github.io/),
    Princeton University
