import sys
import time
import networkx as nx
import numpy as np
import pandas as pd
sys.path.insert(0, "../python")

from motifcluster import motifadjacency as mcmo
from motifcluster import sampling as mcsa


# functions
#######################################################################

def performance_trial(ns, k, motifs, method, nreps, graph_type):

  results = pd.DataFrame(columns = ["n", "k", "motif", "method", "time"])

  for n in reversed(sorted(ns)):
    for motif in motifs:
      for rep in range(nreps):

        # graph parameters
        block_sizes = [n]
        connection_matrix = np.array([k/n]).reshape((1, 1))

        # sample graph
        if graph_type == "erdos_renyi":
          adj_mat = mcsa.sample_dsbm(block_sizes, connection_matrix)
        elif graph_type == "barabasi_albert":
          sample_graph = nx.barabasi_albert_graph(n, k)
          adj_mat = nx.to_scipy_sparse_array(sample_graph)

        if method == "dense":
          adj_mat = adj_mat.toarray()

        # time mam construction
        t0 = time.time()
        mcmo.build_motif_adjacency_matrix(adj_mat, motif, "func", "mean", method)
        t1 = time.time()
        dt = t1 - t0

        # add results
        results_dict = {"n":n, "k":k, "motif":motif, "method": method,
                        "rep":rep, "time":dt}
        results_frame = pd.Series(results_dict).to_frame().T
        results = pd.concat([results, results_frame], ignore_index=True)

        print("n = ",n, ", k = ",k, ", ", motif, ", method = ", method,
              ", rep = ", rep, ", graph type = ", graph_type,
              ", time = {:.3f}".format(dt), sep="")

  # save results
  results.to_csv("results/python_k" + str(k) + "_" + method +
                 "_" + graph_type + ".csv", index=False)

  return None



# script
#######################################################################

motifs = ['M1','M8','M11']
nreps = 10

ns= [101, 200, 500, 1000]
performance_trial(ns, 100, motifs, "dense", nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, "dense", nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, "dense", nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, "dense", nreps, "erdos_renyi")

ns= [101, 200, 500, 1000]
performance_trial(ns, 100, motifs, "sparse", nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, "sparse", nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, "sparse", nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, "sparse", nreps, "erdos_renyi")
