"""
Assorted utility functions for the motifcluster module
are in `motifcluster.utils`.
"""

import random
import numpy as np
from numpy import random as rd
import networkx as nx
from scipy import sparse

def _a_b_one(a_mat, b_mat):

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
  sparse matrix
    The sparse square matrix `a * (b @ one_mat)`.
  """

  if not sparse.issparse(a_mat):
    a_mat = sparse.csr_matrix(a_mat)

  if not sparse.issparse(b_mat):
    b_mat = sparse.csr_matrix(b_mat)

  n = a_mat.shape[0]
  ones_vec = np.ones(n)
  ans = (a_mat.T.multiply(b_mat @ ones_vec)).T

  return ans


def _a_one_b(a_mat, b_mat):

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
  sparse matrix
    The sparse square matrix `a * (one_mat @ b)`.
  """

  if not sparse.issparse(a_mat):
    a_mat = sparse.csr_matrix(a_mat)

  if not sparse.issparse(b_mat):
    b_mat = sparse.csr_matrix(b_mat)

  n = a_mat.shape[0]
  ones_vec = np.ones(n)
  ans = (a_mat.multiply(ones_vec @ b_mat))

  return ans


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
  sparse_mat : sparse matrix
    A sparse-form copy of `some_mat` with its
    diagonal entries set to zero.
  """


  if sparse.issparse(some_mat):
    I = sparse.identity(some_mat.shape[0])
    ans = some_mat - I.multiply(some_mat)

  else:
    I = np.identity(some_mat.shape[0])
    ans = some_mat - some_mat * I

  return ans


def get_largest_component(adj_mat, gr_method):

  """
  Get largest connected component.

  Get the indices of the vertices in the largest connected
  component of a graph from its adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    An adjacency matrix of a graph.
  gr_method : str
    Format to use before building the graph.
    One of `"sparse"` or `"dense"`.

  Returns
  -------
  verts_to_keep : list
    A list of indices corresponding to the vertices in the largest
    connected component.

  Examples
  --------
  >>> adj_mat = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0]).reshape((3, 3))
  >>> get_largest_component(adj_mat)
  """

  if gr_method == "sparse":

    if not sparse.issparse(adj_mat):
      adj_mat = sparse.csr_matrix(adj_mat)

    gr = nx.from_scipy_sparse_array(adj_mat > 0)

  else:

    if isinstance(adj_mat, np.ndarray): # pylint: disable=else-if-used
      gr = nx.from_numpy_array(1 * np.array(adj_mat > 0))

    else:
      gr = nx.from_numpy_array(1 * (adj_mat > 0).toarray())

  verts_to_keep = max(nx.connected_components(gr), key=len)
  verts_to_keep = sorted(verts_to_keep)

  return verts_to_keep


def get_motif_names():

  """
  Get common motif names.

  Get the names of some common motifs as strings.

  Returns
  -------
  motif_names : list
    A list of names (strings) of common motifs.
  """

  motif_names = ["Ms", "Md"]

  for i in range(1, 14):
    motif_name = "M" + str(i)
    motif_names = motif_names + [motif_name]

  motif_names = motif_names + ["Mcoll"] + ["Mexpa"]

  return motif_names


def _random_sparse_matrix(m, n, p, sample_weight_type="constant", w=1):

  """
  Build a random sparse matrix.

  Build a sparse matrix of size `m * n` with non-zero probability `p`.
  Edge weights can be unweighted, constant-weighted or
  Poisson-weighted.

  Parameters
  ----------
  m : int
    Size of first dimension of matrix.
  n : int
    Size of second dimension of matrix.
  p : float
    Probability that each entry is non-zero (before weighting).
  sample_weight_type : str
    Type of weighting scheme.
  w : float
    Weight parameter.

  Returns
  -------
  sparse matrix
    A random sparse matrix.
  """

  mn = m * n

  # number of nonzero entries
  k = rd.binomial(mn, p)

  # indices of nonzero entries
  zs = k * [0]
  inds = random.sample(range(mn), k)

  # values to go in matrix
  if sample_weight_type == "constant":
    vals = np.full(k, w)

  elif sample_weight_type == "poisson":
    vals = rd.poisson(w, size=k)

  else:
    vals = np.ones(k)

  # create small matrix
  if mn <= 1e4:
    ans = np.zeros(mn)
    ans[np.array(inds, dtype=int)] = np.array(vals, dtype=int)
    ans = ans.reshape((m, n))

  # create large matrix
  else:
    ans = sparse.csc_matrix((vals, (inds, zs)), shape=(mn, 1))
    ans = ans.reshape((m, n))

  return ans
