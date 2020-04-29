import sys
import time
sys.path.append("..")

from motifcluster import clustering as mccl
from motifcluster import motifadjacency as mcmo
from motifcluster import sampling as mcsa
from motifcluster import utils as mcut

import numpy as np

def test_perf():

  for n in [100, 200, 500, 1000, 2000]:
    for p in [10/n, 100/n]:
      for motif_name in ["Ms", "Md", "M1", "M9", "M11"]:

        t0 = time.time()

        block_sizes = [n]
        connection_matrix = np.array([p]).reshape((1, 1))
        weight_matrix = np.array([10]).reshape((1, 1))

        adj_mat = mcsa.sample_dsbm(block_sizes, connection_matrix,
                                   weight_matrix)

        mcmo.build_motif_adjacency_matrix(adj_mat, motif_name)

        #print(n, p, motif_name)
        #print("{:.2f}".format(time.time() - t0), "\n")

  return None
