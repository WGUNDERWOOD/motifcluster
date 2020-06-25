import sys
import time
sys.path.insert(0, "../python")

from motifcluster import motifadjacency as mcmo
from motifcluster import sampling as mcsa

import numpy as np
import pandas as pd

# functions
#######################################################################

def performance_trial(ns, k, motifs, method, nreps):

  results = pd.DataFrame(columns = ["n", "k", "motif", "method", "time"])

  for n in reversed(sorted(ns)):
    for motif in motifs:
      for rep in range(nreps):

        # sample graph
        block_sizes = [n]
        connection_matrix = np.array([k/n]).reshape((1, 1))
        adj_mat = mcsa.sample_dsbm(block_sizes, connection_matrix)
        if method == "dense":
          adj_mat = adj_mat.toarray()

        # time mam construction
        t0 = time.time()
        mcmo.build_motif_adjacency_matrix(adj_mat, motif, "func", "mean", method)
        t1 = time.time()
        dt = t1 - t0

        # add results
        results_dict = {"n":n, "k":k, "motif":motif, "method": method, "rep":rep, "time":dt}
        results = results.append(results_dict, ignore_index=True)

        print("n = ",n, ", k = ",k, ", ", motif, ", method = ", method,
              ", rep = ", rep, ", time = {:.3f}".format(dt), sep="")

  # save results
  results.to_csv("results/python_k" + str(k) + "_" + method + ".csv", index=False)

  return None



# script
#######################################################################

motifs = ['M1','M8','M9','M11']
nreps = 100


#ns= [100, 200, 500, 1000, 2000]
ns= [100, 200]
performance_trial(ns, 10, motifs, "dense", nreps)

#performance_trial(ns, 10, motifs, "sparse", nreps)

#performance_trial(ns, 100, motifs, "dense", nreps)

#performance_trial(ns, 100, motifs, "sparse", nreps)
