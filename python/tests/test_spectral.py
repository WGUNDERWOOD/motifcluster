from motifcluster import spectral as mcsp
import numpy as np

# get_first_eigs

def test_get_first_eigs_dense():

  G = np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4))

  vals = [-9, 18, 27]
  vects = np.array([-2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0]).reshape((3, 4)).transpose() / 3

  spect = mcsp.get_first_eigs(G, 3)

  assert np.allclose(spect["vals"], vals, atol = 1e-6)
  assert np.allclose(spect["vects"], vects, atol = 1e-6)

#test_that("get_first_eigs returns correct values on sparse matrix", {

  #G = drop0(matrix(c(7, -4, 14, 0, -4, 19, 10, 0,
                      #14, 10, 10, 0, 0, 0, 0, 100), nrow = 4))
  #vals = c(-9, 18, 27)
  #vects = matrix(c(-2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0) / 3, nrow = 4)

  #spect = mcsp.get_first_eigs(G, 3)

  #expect_equal(spect$vals, vals)
  #expect_equal(spect$vects, vects)
#})

# build_laplacian

#test_that("build_laplacian returns correct matrices on dense matrix", {

  #G = matrix(c(0:8), nrow = 3)
  #G = drop0(G + t(G))

  #degs_mat = diag(c(12, 24, 36))
  #comb_lap = degs_mat - G
  #rw_lap = drop0(solve(degs_mat) %*% (degs_mat - G))

  #expect_equal(build_laplacian(G, type_lap = "comb"), comb_lap)
  #expect_equal(build_laplacian(G, type_lap = "rw"), rw_lap)
#})

#test_that("build_laplacian returns correct matrices on sparse matrix", {

  #G = drop0(matrix(c(0:8), nrow = 3))
  #G = G + t(G)

  #degs_mat = diag(c(12, 24, 36))
  #comb_lap = degs_mat - G
  #rw_lap = drop0(solve(degs_mat) %*% (degs_mat - G))

  #expect_equal(build_laplacian(G, type_lap = "comb"), comb_lap)
  #expect_equal(build_laplacian(G, type_lap = "rw"), rw_lap)
#})

#test_that("build_laplacian gives correct error if row sums are zero", {

  #G = drop0(matrix(c(0, 1, 0, 2)))
  #expect_error(build_laplacian(G, type_lap = "rw"),
               #"row sums of adj_mat must be non-zero")
#})
