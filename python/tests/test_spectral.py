from motifcluster import spectral as mcsp

import random
from scipy import sparse
import numpy as np
from numpy import linalg
import pytest

# get_first_eigs

def test_get_first_eigs_dense():

  G = np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4))

  ans_vals = [-9, 18, 27]
  ans_vects = np.array([-2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0]).reshape((3, 4)).transpose() / 3

  spect = mcsp._get_first_eigs(G, 3)
  vects = spect["vects"]
  vals = spect["vals"]

  for i in range(len(vals)):
    if np.sign(vects[0, i]) != np.sign(ans_vects[0, i]):
      vects[:, i] = -vects[:, i]

  assert np.allclose(vals, ans_vals, atol = 1e-6)
  assert np.allclose(vects, ans_vects, atol = 1e-6)


def test_get_first_eigs_sparse():

  G = sparse.csr_matrix(np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4)), dtype="float")

  ans_vals = [-9, 18, 27]
  ans_vects = np.array([2, 1, -2, 0, -2, 2, -1, 0, 1, 2, 2, 0]).reshape((3, 4)).transpose() / 3

  spect = mcsp._get_first_eigs(G, 3)
  vects = spect["vects"]
  vals = spect["vals"]

  for i in range(len(vals)):
    if np.sign(vects[0, i]) != np.sign(ans_vects[0, i]):
      vects[:, i] = -vects[:, i]

  assert np.allclose(vals, ans_vals, atol = 1e-6)
  assert np.allclose(vects, ans_vects, atol = 1e-6)


def test_get_few_first_eigs_dense():

  G = np.array([7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100]).reshape((4, 4))

  ans_vals = [-9, 18]

  ans_vects = np.array([2, 1, -2, 0, -2, 2, -1, 0]).reshape((2, 4)).transpose() / 3

  spect = mcsp._get_first_eigs(G, 2)
  vects = spect["vects"]
  vals = spect["vals"]

  for i in range(len(vals)):
    if np.sign(vects[0, i]) != np.sign(ans_vects[0, i]):
      vects[:, i] = -vects[:, i]

  assert np.allclose(vals, ans_vals, atol = 1e-6)
  assert np.allclose(vects, ans_vects, atol = 1e-6)


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


# run_laplace_embedding

def test_run_laplace_embedding_dense():

  np.random.seed(9235)

  G = np.array(range(9)).reshape((3, 3))
  G = G + G.transpose()

  ans_vals_comb = [0, 17.07]
  ans_vects_comb = np.array([
    0.577,  0.789,
    0.577, -0.577,
    0.577, -0.211
  ]).reshape((3, 2))

  ans_vals_rw = [0, 1]
  ans_vects_rw = np.array([
    0.577,  0.408,
    0.577, -0.816,
    0.577,  0.408
  ]).reshape((3, 2))

  spectrum_comb = mcsp.run_laplace_embedding(G, 2, "comb")
  spectrum_rw = mcsp.run_laplace_embedding(G, 2, "rw")

  vals_comb = spectrum_comb["vals"]
  vects_comb = spectrum_comb["vects"]
  vals_rw = spectrum_rw["vals"]
  vects_rw = spectrum_rw["vects"]

  for i in range(len(vals_comb)):
    if np.sign(vects_comb[0, i]) != np.sign(ans_vects_comb[0, i]):
      vects_comb[:, i] = -vects_comb[:, i]
    if np.sign(vects_rw[0, i]) != np.sign(ans_vects_rw[0, i]):
      vects_rw[:, i] = -vects_rw[:, i]

  assert np.allclose(vals_comb, ans_vals_comb, atol=0.01)
  assert np.allclose(vects_comb, ans_vects_comb, atol=0.01)
  assert np.allclose(vals_rw, ans_vals_rw, atol=0.01)
  assert np.allclose(vects_rw, ans_vects_rw, atol=0.01)


def test_run_laplace_embedding_sparse():

  np.random.seed(9235)

  G = np.array(range(9)).reshape((3, 3))
  G = sparse.csr_matrix(G + G.transpose())

  ans_vals_comb = [0, 17.07]
  ans_vects_comb = np.array([
    0.577,  0.789,
    0.577, -0.577,
    0.577, -0.211
  ]).reshape((3, 2))

  ans_vals_rw = [0, 1]
  ans_vects_rw = np.array([
    0.577,  0.408,
    0.577, -0.816,
    0.577,  0.408
  ]).reshape((3, 2))

  spectrum_comb = mcsp.run_laplace_embedding(G, 2, "comb")
  spectrum_rw = mcsp.run_laplace_embedding(G, 2, "rw")

  vals_comb = spectrum_comb["vals"]
  vects_comb = spectrum_comb["vects"]
  vals_rw = spectrum_rw["vals"]
  vects_rw = spectrum_rw["vects"]

  for i in range(len(vals_comb)):
    if np.sign(vects_comb[0, i]) != np.sign(ans_vects_comb[0, i]):
      vects_comb[:, i] = -vects_comb[:, i]
    if np.sign(vects_rw[0, i]) != np.sign(ans_vects_rw[0, i]):
      vects_rw[:, i] = -vects_rw[:, i]

  assert np.allclose(vals_comb, ans_vals_comb, atol=0.01)
  assert np.allclose(vects_comb, ans_vects_comb, atol=0.01)
  assert np.allclose(vals_rw, ans_vals_rw, atol=0.01)
  assert np.allclose(vects_rw, ans_vects_rw, atol=0.01)


# run_motif_embedding

def test_run_mot_embedding_dense_restrict():

  np.random.seed(9235)

  adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4))

  # answers
  ans_adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4))

  ans_motif_adj_mat = np.array([
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4))

  ans_comps = [0, 1, 2]

  ans_adj_mat_comps = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3))

  ans_motif_adj_mat_comps = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ]).reshape((3, 3))

  ans_vals = [0, 1.354]

  ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
  ]).reshape((3, 2))

  # run motif embedding
  emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                      "rw", restrict=True)

  # flip eigenvector signs if necessary
  for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
      emb_list["vects"][:, i] = -emb_list["vects"][:, i]

  assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
  assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
  assert np.allclose(ans_comps, emb_list["comps"])
  assert np.allclose(ans_adj_mat_comps, emb_list["adj_mat_comps"].toarray())
  assert np.allclose(ans_motif_adj_mat_comps, emb_list["motif_adj_mat_comps"].toarray())
  assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
  assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


def test_run_mot_embedding_sparse_restrict():

  np.random.seed(9235)

  adj_mat = sparse.csr_matrix(np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4)))

  # answers
  ans_adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4))

  ans_motif_adj_mat = np.array([
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
  ]).reshape((4, 4))

  ans_comps = [0, 1, 2]

  ans_adj_mat_comps = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3))

  ans_motif_adj_mat_comps = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ]).reshape((3, 3))

  ans_vals = [0, 1.354]

  ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
  ]).reshape((3, 2))

  # run motif embedding
  emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                      "rw", restrict=True)

  # flip eigenvector signs if necessary
  for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
      emb_list["vects"][:, i] = -emb_list["vects"][:, i]

  assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
  assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
  assert np.allclose(ans_comps, emb_list["comps"])
  assert np.allclose(ans_adj_mat_comps, emb_list["adj_mat_comps"].toarray())
  assert np.allclose(ans_motif_adj_mat_comps, emb_list["motif_adj_mat_comps"].toarray())
  assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
  assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


def test_run_mot_embedding_dense_no_restrict():

  np.random.seed(9235)

  adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3))

  # answers
  ans_adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3))

  ans_motif_adj_mat = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ]).reshape((3, 3))

  ans_vals = [0, 1.354]

  ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
  ]).reshape((3, 2))

  # run motif embedding
  emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                      "rw", restrict=False)

  # flip eigenvector signs if necessary
  for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
      emb_list["vects"][:, i] = -emb_list["vects"][:, i]

  assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
  assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
  assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
  assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


def test_run_mot_embedding_sparse_no_restrict():

  np.random.seed(9235)

  adj_mat = sparse.csr_matrix(np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3)))

  # answers
  ans_adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ]).reshape((3, 3))

  ans_motif_adj_mat = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ]).reshape((3, 3))

  ans_vals = [0, 1.354]

  ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
  ]).reshape((3, 2))

  # run motif embedding
  emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                      "rw", restrict=False)

  # flip eigenvector signs if necessary
  for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
      emb_list["vects"][:, i] = -emb_list["vects"][:, i]

  assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
  assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
  assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
  assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)
