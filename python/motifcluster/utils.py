"""
Assorted utility functions for the motifcluster module
are in `motifcluster.utils`.
"""

import numpy as np
from numpy import random as rd
import networkx as nx
from scipy import sparse

def _a_b_one(a, b):

  """
  Compute a right-multiplication with the ones matrix.

  Compute `a * (b @ one_mat)` where `a`, `b`,
  `ones_mat` are square matrices of the same size,
  and `ones_mat` contains all entries equal to one.
  The product `*` is an entry-wise (Hadamard) product,
  while `@` represents matrix multiplication.
  This method is more efficient than the naive approach
  when `a` or `b` are sparse.

  Parameters
  ----------
  a, b : matrix
    Square matrices of the same size.

  Returns
  -------
  scipy.sparse.csr_matrix
    The sparse square matrix `a * (b @ one_mat)`.
  """

  a_sparse = sparse.csr_matrix(a)
  b_sparse = sparse.csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.T.multiply(b_sparse @ ones_vec)).T

  return sparse.csr_matrix(ans)


def _a_one_b(a, b):

  """
  Compute a left-multiplication with the ones matrix.

  Compute `a * (one_mat @ b)` where `a`, `b`,
  `ones_mat` are square matrices of the same size,
  and `ones_mat` contains all entries equal to one.
  The product `*` is an entry-wise (Hadamard) product,
  while `@` represents matrix multiplication.
  This method is more efficient than the naive approach
  when `a` or `b` are sparse.

  Parameters
  ----------
  a, b : matrix
    Square matrices of the same size.

  Returns
  -------
  scipy.sparse.csr_matrix
    The sparse square matrix `a * (one_mat @ b)`.
  """

  a_sparse = sparse.csr_matrix(a)
  b_sparse = sparse.csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.multiply(ones_vec @ b))

  return sparse.csr_matrix(ans)


def _drop0_killdiag(some_mat):

  """
  Set diagonal entries to zero and sparsify.

  Set the diagonal entries of a matrix to zero
  and convert it to sparse form.

  Parameters
  ----------
  some_mat : matrix
    A square matrix.

  Returns
  -------
  scipy.sparse.csr_matrix
    A sparse-form copy of `some_matrix` with its
    diagonal entries set to zero.
  """

  ans = sparse.csr_matrix(some_mat)
  I = sparse.identity(ans.shape[0])
  ans = ans - I.multiply(ans)

  return(sparse.csr_matrix(ans))


def get_largest_component(adj_mat):

  """
  Get largest connected component.

  Get the indices of the vertices in the largest connected
  component of a graph from its adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    An adjacency matrix of a graph.

  Returns
  -------
  list
    A list of indices corresponding to the vertices in the largest
    connected component.

  Examples
  --------
  >>> adj_mat = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0]).reshape((3, 3))
  >>> get_largest_component(adj_mat)
  """

  adj_mat_sparse = sparse.csr_matrix(adj_mat)
  gr = nx.from_scipy_sparse_matrix(adj_mat_sparse > 0)
  verts_to_keep = max(nx.connected_components(gr), key=len)
  verts_to_keep = sorted(verts_to_keep)

  return verts_to_keep


def get_motif_names():

  """
  Get common motif names.

  Get the names of some common motifs as strings.

  Returns
  -------
  list
    A list of names (strings) of common motifs.
  """

  motif_names = ["Ms", "Md"]

  for i in range(1, 14):
    motif_name = "M" + str(i)
    motif_names = [motif_names, motif_name]

  motif_names = [motif_names, "Mcoll", "Mexpa"]

  return(motif_names)


def random_sparse_matrix(m, n, p):

  # TODO doc and test
  mn = m * n

  # number of nonzero entries
  k = rd.binomial(mn, p)
  vals = np.ones(k)

  # indices of nonzero entries
  inds = rd.choice(mn, size = k, replace = False)
  zs = k * [0]

  ans = sparse.csr_matrix((vals, (inds, zs)), shape = (mn, 1))
  ans = ans.reshape((m, n))

  return(ans)
