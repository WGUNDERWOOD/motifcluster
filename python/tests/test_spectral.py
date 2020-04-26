from motifcluster import spectral as mcsp
from scipy import sparse
import numpy as np
from numpy import linalg
import pytest

# get_first_eigs

def test_get_first_eigs_dense():

  G = np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4))

  vals = [-9, 18, 27]
  vects = np.array([-2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0]).reshape((3, 4)).transpose() / 3

  spect = mcsp.get_first_eigs(G, 3)

  assert np.allclose(spect["vals"], vals, atol = 1e-6)
  assert np.allclose(spect["vects"], vects, atol = 1e-6)


def test_get_first_eigs_sparse():

  G = sparse.csr_matrix(np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4)))

  vals = [-9, 18, 27]
  vects = np.array([2, 1, -2, 0, -2, 2, -1, 0, 1, 2, 2, 0]).reshape((3, 4)).transpose() / 3

  spect = mcsp.get_first_eigs(G, 3)
  print(vects)
  print(spect["vects"])

  assert np.allclose(spect["vals"], vals, atol = 1e-6)
  assert np.allclose(spect["vects"], vects, atol = 1e-6)


# build_laplacian

def test_build_laplacian_dense():

  G = np.array(range(9)).reshape((3, 3))
  G = G + G.transpose()

  degs_mat = np.diag([12, 24, 36])
  comb_lap = degs_mat - G
  rw_lap = linalg.inv(degs_mat) @ (degs_mat - G)

  assert np.allclose(mcsp.build_laplacian(G, type_lap = "comb").toarray(), comb_lap)
  assert np.allclose(mcsp.build_laplacian(G, type_lap = "rw").toarray(), rw_lap)


def test_build_laplacian_sparse():

  G = np.array(range(9)).reshape((3, 3))
  G = sparse.csr_matrix(G + G.transpose())

  degs_mat = np.diag([12, 24, 36])
  comb_lap = degs_mat - G
  rw_lap = linalg.inv(degs_mat) @ (degs_mat - G)

  assert np.allclose(mcsp.build_laplacian(G, type_lap = "comb").toarray(), comb_lap)
  assert np.allclose(mcsp.build_laplacian(G, type_lap = "rw").toarray(), rw_lap)


def test_build_laplacian_row_sum_error():

  with pytest.raises(AssertionError):
    G = np.array([0, 1, 0, 2]).reshape((2,2))
    mcsp.build_laplacian(G, type_lap = "rw")
