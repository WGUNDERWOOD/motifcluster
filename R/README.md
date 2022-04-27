# motifcluster <img src="https://github.com/WGUNDERWOOD/motifcluster/raw/develop/sticker/hex_sticker_small.png" alt="motifcluster sticker" align="right" width=120/>

An R package for motif-based spectral clustering of weighted directed networks.

## Introduction

The **motifcluster** package provides
implementations of motif-based spectral clustering
of weighted directed networks in R.
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
install.packages("motifcluster")
```

## Dependencies

- ClusterR
- igraph
- Matrix
- RSpectra

## Documentation

Documentation for the **motifcluster** package
is available in the
[doc](https://github.com/WGUNDERWOOD/motifcluster/tree/master/R/doc)
directory.

## Vignette

An instructional vignette for the **motifcluster** package
is available in the
[vignettes](https://github.com/WGUNDERWOOD/motifcluster/tree/master/R/vignettes)
directory.

## Author

- [William George Underwood](https://wgunderwood.github.io/),
    Princeton University
