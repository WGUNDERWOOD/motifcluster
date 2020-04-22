from motifcluster import sampling as mcsa

import numpy as np
from numpy import random as rd
from scipy import sparse

def test_sample_dsbm_unweighted():

  rd.seed(seed = 9387)

  sample_weight_type = "unweighted"
  block_sizes = [2, 3]
  connection_matrix = np.array([0.4, 0.5, 0.6, 0.7]).reshape((2, 2))
  n_reps = 200
  weight_matrix = None

  G = mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)

  assert (G.toarray() * (1 - G.toarray()) == 0).all

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

  print(G)
  print(ans)

  assert np.allclose(G, ans, atol = 0.05)
