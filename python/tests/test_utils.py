# import motifcluster
from motifcluster import utils as mcut

# import other dependencies
import numpy as np
from scipy.sparse import csr_matrix

# tests

def test_a_b_one():

  a_dense = np.array(range(-4, 5)).reshape((3, 3))
  b_dense = np.array(range(-1, 8)).reshape((3, 3))
  a_sparse = csr_matrix(a_dense)
  b_sparse = csr_matrix(b_dense)

  ans = np.array([0, 0, 0, -9, 0, 9, 36, 54, 72]).reshape((3, 3))
  ans = csr_matrix(ans)

  assert np.allclose(
    mcut.a_b_one(a_dense, b_dense).toarray(),
    ans.toarray()
  )

  assert np.allclose(
    mcut.a_b_one(a_sparse, b_sparse).toarray(),
    ans.toarray()
  )


def test_a_one_b():

  a_dense = np.array(range(-4, 5)).reshape((3, 3))
  b_dense = np.array(range(-1, 8)).reshape((3, 3))
  a_sparse = csr_matrix(a_dense)
  b_sparse = csr_matrix(b_dense)

  ans = np.array([-24, -27, -24, -6, 0, 12, 12, 27, 48]).reshape((3, 3))
  ans = csr_matrix(ans)

  assert np.allclose(
    mcut.a_one_b(a_dense, b_dense).toarray(),
    ans.toarray()
  )

  assert np.allclose(
    mcut.a_one_b(a_sparse, b_sparse).toarray(),
    ans.toarray()
  )
