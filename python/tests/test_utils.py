from motifcluster import utils as mcut

import numpy as np
from scipy import sparse

def test_a_b_one():

  a_dense = np.array(range(-4, 5)).reshape((3, 3))
  b_dense = np.array(range(-1, 8)).reshape((3, 3))
  a_sparse = sparse.csr_matrix(a_dense)
  b_sparse = sparse.csr_matrix(b_dense)

  ans = np.array([0, 0, 0, -9, 0, 9, 36, 54, 72]).reshape((3, 3))

  assert np.allclose(
    mcut._a_b_one(a_dense, b_dense).toarray(),
    ans
  )

  assert np.allclose(
    mcut._a_b_one(a_sparse, b_sparse).toarray(),
    ans
  )


def test_a_one_b():

  a_dense = np.array(range(-4, 5)).reshape((3, 3))
  b_dense = np.array(range(-1, 8)).reshape((3, 3))
  a_sparse = sparse.csr_matrix(a_dense)
  b_sparse = sparse.csr_matrix(b_dense)

  ans = np.array([-24, -27, -24, -6, 0, 12, 12, 27, 48]).reshape((3, 3))

  assert np.allclose(
    mcut._a_one_b(a_dense, b_dense).toarray(),
    ans
  )

  assert np.allclose(
    mcut._a_one_b(a_sparse, b_sparse).toarray(),
    ans
  )


def test_drop0_kill_diag():

  adj_mat_dense = np.array(range(-1, 8)).reshape((3,3))
  adj_mat_sparse = sparse.csr_matrix(adj_mat_dense)

  ans = np.array([0, 0, 1, 2, 0, 4, 5, 6, 0]).reshape((3,3))

  assert np.allclose(
    mcut._drop0_killdiag(adj_mat_dense),
    ans
  )

  assert np.allclose(
    mcut._drop0_killdiag(adj_mat_sparse).toarray(),
    ans
  )


def test_get_largest_component():

  adj_mat_dense = np.array([
    0, 0, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 0, 2,
    3, 0, 0, 0, 0,
    0, 4, 0, 0, 0
  ]).reshape((5, 5))

  adj_mat_sparse = sparse.csr_matrix(adj_mat_dense)
  ans = [1, 2, 4]

  assert mcut.get_largest_component(adj_mat_dense, "dense") == ans
  assert mcut.get_largest_component(adj_mat_dense, "sparse") == ans
  assert mcut.get_largest_component(adj_mat_sparse, "dense") == ans
  assert mcut.get_largest_component(adj_mat_sparse, "sparse") == ans
