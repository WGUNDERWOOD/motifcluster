
from motifcluster import sampling as mcsa

import numpy as np
from numpy import random as rd
from scipy import sparse
import random

#TODO remove
def test_dsbm_speed():

  sample_weight_type = "unweighted"
  block_sizes = [int(1e3), int(1e3)]
  connection_matrix = np.array(4 * [0.0001]).reshape((2, 2))
  weight_matrix = np.array([4, 4, 4, 4]).reshape((2, 2))
  sample_weight_type = "poisson"

  G = mcsa.sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)


  assert 0 == 0
  return
