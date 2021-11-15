..  image:: https://github.com/WGUNDERWOOD/motifcluster/raw/develop/sticker/hex_sticker_small.png
    :width: 130
    :alt: motifcluster sticker
    :align: right

motifcluster
============


A Python package for motif-based spectral clustering of weighted directed networks.

Introduction
------------
The **motifcluster** package provides
implementations of motif-based spectral clustering
of weighted directed networks in Python.
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
`arxiv:2004.01293 <https://arxiv.org/abs/2004.01293>`_.

Installation
------------

From PyPI:

.. code-block:: bash

  pip install motifcluster

With conda:

.. code-block:: bash

  conda install -c conda-forge motifcluster

Dependencies
------------
- Networkx
- Numpy
- Scipy
- Scikit-learn

Documentation
-------------
Documentation for the **motifcluster** package
is available on
`Read the Docs <https://motifcluster.readthedocs.io/en/latest/>`_.

Tutorial
--------
A tutorial for the **motifcluster** package
is available on Github in the
`tutorial directory <https://github.com/WGUNDERWOOD/motifcluster/blob/develop/python/tutorial/motifcluster_tutorial.pdf>`_.

Authors
-------
- `William George Underwood <https://wgunderwood.github.io/>`_,
  Princeton University (maintainer)
- `Andrew Elliott <https://www.turing.ac.uk/people/researchers/andrew-elliott>`_,
  The Alan Turing Institute

Links
-----
- Source code repository on
  `GitHub <https://github.com/WGUNDERWOOD/motifcluster>`_
- Package index page on
  `PyPI <https://pypi.org/project/motifcluster/>`_
- Documentation on
  `Read the Docs <https://motifcluster.readthedocs.io/en/latest/>`_
