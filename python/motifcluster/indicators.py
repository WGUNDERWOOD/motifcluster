"""
Functions for building adjacency and indicator matrices
are in `motifcluster.indicators`.
"""

import motifcluster.utils as mcut

import numpy as np
from scipy import sparse

def _build_G(adj_mat):

  """
  Build sparse adjacency matrix.

  Build the sparse adjacency matrix `G` from a graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    The adjacency matrix `G` in sparse form.
  """

  G = sparse.csr_matrix(adj_mat)
  return(G)


def _build_J(adj_mat):

  """
  Build directed indicator matrix.

  Build the sparse directed indicator matrix `J`
  from a graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A directed indicator matrix `J` in sparse form.
  """

  G = _build_G(adj_mat)
  J = mcut._drop0_killdiag(G > 0)
  return(J)


def _build_Gs(adj_mat):

  """
  Build single-edge indicator matrix.

  Build the sparse single-edge adjacency matrix `Gs` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A single-edge adjacency matrix `Gs` in sparse form.
  """

  G = _build_G(adj_mat)
  J = _build_J(adj_mat)
  Gs = mcut._drop0_killdiag(G - G.multiply(J.T))
  return(Gs)


def _build_Js(adj_mat):

  """
  Build single-edge indicator matrix.

  Build the sparse single-edge indicator matrix `Js` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A single-edge indicator matrix `Js` in sparse form.
  """

  Gs = _build_Gs(adj_mat)
  Js = mcut._drop0_killdiag(Gs > 0)
  return(Js)


def _build_Gd(adj_mat):

  """
  Build double-edge adjacency matrix.

  Build the sparse double-edge adjacency matrix `Gd` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A double-edge adjacency matrix `Gd` in sparse form.
  """

  J = _build_J(adj_mat)
  G = _build_G(adj_mat)
  Gd = mcut._drop0_killdiag((G + G.T) * J * J.T)
  Gd = mcut._drop0_killdiag((G + G.T).multiply(J).multiply(J.T))
  return(Gd)


def _build_Jd(adj_mat):

  """
  Build double-edge indicator matrix.

  Build the sparse double-edge indicator matrix `Jd` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A double-edge indicator matrix `Jd` in sparse form.
  """

  Gd = _build_Gd(adj_mat)
  Jd = mcut._drop0_killdiag(Gd > 0)
  return(Jd)


def _build_J0(adj_mat):

  """
  Build missing-edge indicator matrix.

  Build the missing-edge indicator matrix `J0` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A missing-edge indicator matrix `J0`.
  """

  G_dense = _build_G(adj_mat).toarray()
  J0 = mcut._drop0_killdiag((G_dense + G_dense.T) == 0)
  return(J0)


def _build_Jn(adj_mat):

  """
  Build vertex-distinct indicator matrix.

  Build the vertex-distinct indicator matrix `Jn` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A vertex-distinct indicator matrix `Jn`.
  """

  n = adj_mat.shape[0]
  Jn = mcut._drop0_killdiag(np.ones((n, n)))
  return(Jn)


def _build_Id(adj_mat):

  """
  Build identity matrix.

  Build the sparse identity matrix `Id` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    An identity matrix `Id` in sparse form.
  """

  n = adj_mat.shape[0]
  Id = sparse.identity(n)
  return(Id)


def _build_Je(adj_mat):

  """
  Build edge-and-diagonal matrix.

  Build the sparse edge-and-diagonal matrix `Ie` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    An edge-and-diagonal matrix `Ie` in sparse form.
  """

  G = _build_G(adj_mat)
  Id = _build_Id(adj_mat)
  Je = sparse.csr_matrix(Id + ((G + G.T) > 0))
  return(Je)


def _build_Gp(adj_mat):

  """
  Build product matrix.

  Build the sparse product matrix `Gp` from a
  graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    The original adjacency matrix.

  Returns
  -------
  sparse matrix
    A product matrix `Gp` in sparse form.
  """

  G = _build_G(adj_mat)
  Gp = mcut._drop0_killdiag(G.multiply(G.T))
  return(Gp)
