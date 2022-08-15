"""
Functions for random sampling of weighted directed networks
are in `motifcluster.sampling`.
"""

import numpy as np
from scipy import sparse

from motifcluster import utils as mcut

def sample_dsbm(block_sizes, connection_matrix,
                weight_matrix=None,
                sample_weight_type="unweighted"):

  """
  Sample a directed stochastic block model (DSBM).

  Sample the (weighted) adjacency matrix of a (weighted) directed stochastic
  block model (DSBM) with specified parameters.

  Parameters
  ----------
  block_sizes : list of int
    A list containing the size of each block of vertices.
  connection_matrix : matrix
    A matrix containing the block-to-block connection probabilities.
  sample_weight_type : str
    The type of weighting scheme.
    One of `"unweighted"`, `"constant"` or `"poisson"`.
  weight_matrix : matrix
    A matrix containing the block-to-block weight
    parameters.
    Unused for `sample_weight_type = "constant"`.
    Defaults to `None`.

  Returns
  -------
  adj_mat : sparse matrix
    A randomly sampled (weighted) adjacency matrix of a DSBM.

  Examples
  --------
  >>> block_sizes = [10, 10]
  >>> connection_matrix = np.array([0.8, 0.1, 0.1, 0.8]).reshape((2, 2))
  >>> weight_matrix = np.array([10, 3, 3, 10]).reshape((2, 2))
  >>> sample_dsbm(block_sizes, connection_matrix, weight_matrix, "poisson")
  """

  # check args
  assert block_sizes == [int(x) for x in block_sizes]
  assert all(x > 0 for x in block_sizes)
  assert len(block_sizes) == connection_matrix.shape[0]
  assert len(block_sizes) == connection_matrix.shape[1]
  assert (connection_matrix >= 0).all()
  assert (connection_matrix <= 1).all()
  assert sample_weight_type in ["unweighted", "constant", "poisson"]

  if sample_weight_type != "unweighted":
    assert weight_matrix is not None
    assert len(block_sizes) == weight_matrix.shape[0]
    assert len(block_sizes) == weight_matrix.shape[1]
    assert (weight_matrix >= 0).all()

  # initialize variables
  k = len(block_sizes)
  block_list = []

  for i in range(k):

    row_list = []

    for j in range(k):

      # block parameters
      ni = block_sizes[i]
      nj = block_sizes[j]
      p = connection_matrix[i, j]

      # generate block
      if sample_weight_type == "unweighted":
        w = 1

      else:
        w = weight_matrix[i, j]

      block = mcut._random_sparse_matrix(ni, nj, p, sample_weight_type, w)
      row_list.append(block)

    block_list.append(row_list)


  # build from dense blocks
  if isinstance(block_list[0][0], np.ndarray):
    adj_mat = np.block(block_list)

  # build from sparse blocks
  else:
    adj_mat = sparse.bmat(block_list)

  adj_mat = sparse.csr_matrix(adj_mat) # pylint: disable=redefined-variable-type
  adj_mat = mcut._drop0_killdiag(adj_mat)

  return adj_mat


def sample_bsbm(source_block_sizes, dest_block_sizes,
                bipartite_connection_matrix,
                bipartite_weight_matrix=None,
                sample_weight_type="unweighted"):

  """
  Sample a bipartite stochastic block model (BSBM).

  Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic
  block model (BSBM) with specified parameters.

  Parameters
  ----------
  source_block_sizes : list of int
    A list containing the size of each block of source vertices.
  dest_block_sizes : list of int
    A list containing the size of each block of destination vertices.
  bipartite_connection_matrix : matrix
    A matrix containing the source block to destination block
    connection probabilities.
  bipartite_weight_matrix : matrix
    A matrix containing the source block to destination block weight
    parameters. Unused for `sample_weight_type = "constant"`.
    Defaults to `None`.
  sample_weight_type : str
    The type of weighting scheme.
    One of `"unweighted"`, `"constant"` or `"poisson"`.

  Returns
  -------
  adj_mat : sparse matrix
    A randomly sampled (weighted) adjacency matrix of a BSBM.

  Examples
  --------
  >>> source_block_sizes = [10, 10]
  >>> dest_block_sizes = [10, 10, 10]
  >>> bipartite_connection_matrix = np.array([0.8, 0.5, 0.1, 0.1, 0.5, 0.8]).reshape((2, 3))
  >>> bipartite_weight_matrix = np.array([20, 10, 2, 2, 10, 20]).reshape((2, 3))
  >>> sample_bsbm(block_sizes, bipartite_connection_matrix,
  ...   bipartite_weight_matrix, "poisson")
  """

  # check args
  assert source_block_sizes == [int(x) for x in source_block_sizes]
  assert dest_block_sizes == [int(x) for x in dest_block_sizes]
  assert all(x > 0 for x in source_block_sizes)
  assert all(x > 0 for x in dest_block_sizes)
  assert len(source_block_sizes) == bipartite_connection_matrix.shape[0]
  assert len(dest_block_sizes) == bipartite_connection_matrix.shape[1]
  assert (bipartite_connection_matrix >= 0).all()
  assert (bipartite_connection_matrix <= 1).all()
  assert sample_weight_type in ["unweighted", "constant", "poisson"]

  if sample_weight_type != "unweighted":
    assert bipartite_weight_matrix is not None
    assert len(source_block_sizes) == bipartite_weight_matrix.shape[0]
    assert len(dest_block_sizes) == bipartite_weight_matrix.shape[1]
    assert (bipartite_weight_matrix >= 0).all()

  # initialize parameters
  ks = len(source_block_sizes)
  kd = len(dest_block_sizes)
  zeros_ss = np.zeros((ks, ks))
  zeros_d = np.zeros((kd, ks + kd))

  # build block sizes vector
  block_sizes = source_block_sizes + dest_block_sizes

  # build connection matrix
  connection_matrix = np.block([[zeros_ss, bipartite_connection_matrix], [zeros_d]])

  # build weight matrix
  if bipartite_weight_matrix is not None:
    weight_matrix = np.block([[zeros_ss, bipartite_weight_matrix], [zeros_d]])

  else:
    weight_matrix = None

  # sample BSBM
  adj_mat = sample_dsbm(block_sizes, connection_matrix,
                        weight_matrix, sample_weight_type)

  return adj_mat


def demonstration_graph():

  """
  Generate a small graph for demonstrations.

  Generate the sparse and dense adjacency matrices of a small weighted
  directed graph, for demonstrating methods and running tests.

  Returns
  -------
  adj_mat_dense : matrix
    the adjacency matrix in dense form.
  adj_mat_sparse : sparse matrix
    the adjacency matrix in sparse form.
  """

  adj_mat_dense = np.array([
    0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,
    2, 0,  3,  0,  6,  8,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0, 10,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,
    4, 5,  0,  0,  0, 14,  0,  0, 18, 19,  0, 0,
    0, 7,  9,  0, 13,  0,  0,  0,  0,  0, 21, 0,
    0, 0, 11, 12,  0, 15,  0, 17,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 16,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0,  0,  0,  0, 24,  0, 0,
    0, 0,  0,  0,  0, 20,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 22,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 23,  0,  0,  0,  0, 0
  ]).reshape((12, 12))

  adj_mat_sparse = sparse.csr_matrix(adj_mat_dense)

  ans = {
    "adj_mat_dense": adj_mat_dense,
    "adj_mat_sparse": adj_mat_sparse,
  }

  return ans
