from motifcluster import sampling as mcsa

import numpy as np
from numpy import random as rd
from scipy import sparse

# sample_dsbm

def test_sample_dsbm_unweighted():

  rd.seed(seed = 9387)

  sample_weight_type = "unweighted"
  block_sizes = [2, 3]
  connection_matrix = np.array([0.4, 0.5, 0.6, 0.7]).reshape((2, 2))
  n_reps = 200
  weight_matrix = None

  G = mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)

  G_vals = [0, 1]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((5, 5))

  for rep in range(n_reps):
    G += mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                          sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
           0, 0.4, 0.5, 0.5, 0.5,
         0.4,   0, 0.5, 0.5, 0.5,
         0.6, 0.6,   0, 0.7, 0.7,
         0.6, 0.6, 0.7,   0, 0.7,
         0.6, 0.6, 0.7, 0.7,   0
       ]).reshape((5, 5))

  assert np.allclose(G, ans, atol = 0.05)

  return


def test_sample_dsbm_constant_weighted():

  rd.seed(seed = 2839)

  sample_weight_type = "constant"
  block_sizes = [2, 3]
  connection_matrix = np.array([0.4, 0.5, 0.6, 0.7]).reshape((2, 2))
  n_reps = 200
  weight_matrix = np.array([20, 30, 40, 50]).reshape((2, 2))

  G = mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)

  G_vals = [0, 20, 30, 40, 50]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((5, 5))

  for rep in range(n_reps):
    G += mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                          sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
                          0,  8, 15, 15, 15,
                          8,  0, 15, 15, 15,
                         24, 24,  0, 35, 35,
                         24, 24, 35,  0, 35,
                         24, 24, 35, 35,  0
       ]).reshape((5, 5))

  assert np.allclose(G, ans, atol = 3)

  return


def test_sample_dsbm_poisson_weighted():

  rd.seed(seed = 2839)

  sample_weight_type = "poisson"
  block_sizes = [2, 3]
  connection_matrix = np.array([0.4, 0.5, 0.6, 0.7]).reshape((2, 2))
  n_reps = 200
  weight_matrix = np.array([20, 30, 40, 50]).reshape((2, 2))

  G = mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)

  assert (G.toarray() == np.floor(G.toarray())).all()
  assert (G.toarray() >= 0).all()

  G = np.zeros((5, 5))

  for rep in range(n_reps):
    G += mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                          sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
                          0,  8, 15, 15, 15,
                          8,  0, 15, 15, 15,
                         24, 24,  0, 35, 35,
                         24, 24, 35,  0, 35,
                         24, 24, 35, 35,  0
       ]).reshape((5, 5))

  assert np.allclose(G, ans, atol = 3)

  return


# sample_bsbm

def test_sample_bsbm_unweighted():

  rd.seed(seed = 9423)

  sample_weight_type = "unweighted"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 200
  bipartite_weight_matrix = None

  G = mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G_vals = [0, 1]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((6, 6))

  for rep in range(n_reps):
    G += mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
                         0, 0, 0, 0.3, 0.4, 0.5,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 0.05)

  return


def test_sample_bsbm_constant_weighted():

  rd.seed(seed = 7482)

  sample_weight_type = "constant"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 200
  bipartite_weight_matrix = np.array([10, 20, 30, 40, 50, 60]).reshape((2, 3))

  G = mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G_vals = [0, 10, 20, 30, 40, 50, 60]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((6, 6))

  for rep in range(n_reps):
    G += mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 3)

  return


def test_sample_bsbm_poisson_weighted():

  rd.seed(seed = 7482)

  sample_weight_type = "poisson"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 200
  bipartite_weight_matrix = np.array([10, 20, 30, 40, 50, 60]).reshape((2, 3))

  G = mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  assert (G.toarray() == np.floor(G.toarray())).all()
  assert (G.toarray() >= 0).all()

  G = np.zeros((6, 6))

  for rep in range(n_reps):
    G += mcsa.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type) / n_reps

  G = np.array(G)

  ans = np.array([
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 3)

  return
