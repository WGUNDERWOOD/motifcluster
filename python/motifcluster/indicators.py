"""
Functions for building adjacency and indicator matrices
are in `motifcluster.indicators`.
"""

import numpy as np
from scipy import sparse

import motifcluster.utils as mcut

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
  G : sparse matrix
    The adjacency matrix in sparse form.
  """

  G = adj_mat

  return G


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
  J : sparse matrix
    A directed indicator matrix in sparse form.
  """

  G = _build_G(adj_mat)

  if sparse.issparse(adj_mat):
    J = mcut._drop0_killdiag(G > 0)
  else:
    J = mcut._drop0_killdiag(1 * (G > 0))

  return J


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
  Gs : sparse matrix
    A single-edge adjacency matrix in sparse form.
  """

  G = _build_G(adj_mat)
  J = _build_J(adj_mat)

  if sparse.issparse(adj_mat):
    Gs = mcut._drop0_killdiag(G - G.multiply(J.T))
  else:
    Gs = mcut._drop0_killdiag(G - G * J.T)

  return Gs


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
  Js : sparse matrix
    A single-edge indicator matrix in sparse form.
  """

  Gs = _build_Gs(adj_mat)

  if sparse.issparse(adj_mat):
    Js = mcut._drop0_killdiag(Gs > 0)
  else:
    Js = mcut._drop0_killdiag(1 * (Gs > 0))

  return Js


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
  Gd : sparse matrix
    A double-edge adjacency matrix in sparse form.
  """

  J = _build_J(adj_mat)
  G = _build_G(adj_mat)

  if sparse.issparse(adj_mat):
    Gd = mcut._drop0_killdiag((G + G.T).multiply(J).multiply(J.T))
  else:
    Gd = mcut._drop0_killdiag((G + G.T) * J * J.T)

  return Gd


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
  Jd : sparse matrix
    A double-edge indicator matrix in sparse form.
  """

  Gd = _build_Gd(adj_mat)

  if sparse.issparse(adj_mat):
    Jd = mcut._drop0_killdiag(Gd > 0)
  else:
    Jd = mcut._drop0_killdiag(1 * (Gd > 0))

  return Jd


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
  J0 : sparse matrix
    A missing-edge indicator matrix.
  """

  G = _build_G(adj_mat)

  if sparse.issparse(adj_mat):
    J0 = sparse.csr_matrix(mcut._drop0_killdiag((G + G.T) == 0))
  else:
    J0 = mcut._drop0_killdiag((G + G.T) == 0)

  return J0


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
  Jn : sparse matrix
    A vertex-distinct indicator matrix.
  """

  n = adj_mat.shape[0]

  if sparse.issparse(adj_mat):
    Jn = sparse.csr_matrix(mcut._drop0_killdiag(np.ones((n, n))))
  else:
    Jn = mcut._drop0_killdiag(np.ones((n, n)))

  return Jn


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
  Id : sparse matrix
    An identity matrix in sparse form.
  """

  n = adj_mat.shape[0]

  if sparse.issparse(adj_mat):
    Id = sparse.identity(n)
  else:
    Id = np.identity(n)

  return Id


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
  Ie : sparse matrix
    An edge-and-diagonal matrix in sparse form.
  """

  G = _build_G(adj_mat)
  Id = _build_Id(adj_mat)
  Je = Id + ((G + G.T) > 0)

  return Je


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
  Gp : sparse matrix
    A product matrix in sparse form.
  """

  G = _build_G(adj_mat)

  if sparse.issparse(adj_mat):
    Gp = mcut._drop0_killdiag(G.multiply(G.T))
  else:
    Gp = mcut._drop0_killdiag(G * G.T)

  return Gp
