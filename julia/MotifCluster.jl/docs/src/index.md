# MotifCluster

![motifcluster logo](https://github.com/WGUNDERWOOD/motifcluster/raw/main/sticker/hex_sticker_small.png)

A Julia package for motif-based spectral clustering of weighted directed networks.

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
[arxiv:2004.01293](https://arxiv.org/abs/2004.01293).

## Installation

From the Julia General registry:

```
] add MotifCluster
```

## Dependencies
- Aqua
- Clustering
- Distributions
- Graphs

## Documentation
Documentation for the **MotifCluster** package is available on 
[the web](https://wgunderwood.github.io/motifcluster/stable/).

## Tutorial
A tutorial for the **MotifCluster** package
is available on Github in the
[tutorial directory](https://github.com/WGUNDERWOOD/motifcluster/blob/main/julia/MotifCluster.jl/tutorial/motifcluster_tutorial.pdf).

## Authors
- [William George Underwood](https://wgunderwood.github.io/), Princeton University (maintainer)
- [Andrew Elliott](https://www.turing.ac.uk/people/researchers/andrew-elliott), The Alan Turing Institute

## Links
- Source code repository on
  [GitHub](https://github.com/WGUNDERWOOD/motifcluster)
- Documentation on
  [the web](https://wgunderwood.github.io/motifcluster/stable/).
