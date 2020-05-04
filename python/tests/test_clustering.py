import numpy as np
import random
from numpy import random as rd
from sklearn.metrics import adjusted_rand_score

from motifcluster import clustering as mccl
from motifcluster import sampling as mcsa
from motifcluster import utils as mcut

def test_cluster_spectrum():

  rd.seed(2352)
  n = 10
  vects_1 = rd.normal(size = (n, 3))
  vects_2 = np.concatenate((rd.normal(size = (n, 1)),
                            rd.normal(loc = 4, size = (n, 2))),
                            axis = 1)

  vects =  np.concatenate((vects_1, vects_2))
  spectrum = {"vects": vects}

  clust_ans = n*[0] + n*[1]
  clust = mccl.cluster_spectrum(spectrum, 2)

  if clust[1] == 1:
    clust = 1 - clust

  assert (clust == clust_ans).all()


def test_run_motif_clustering():

  rd.seed(3957)
  random.seed(3957)

  n = 50
  block_sizes = 3 * [n]

  connection_matrix = np.array([
    0.9, 0.4, 0.4,
    0.4, 0.9, 0.4,
    0.4, 0.4, 0.9
  ]).reshape((3, 3))

  weight_matrix = np.array([
    9, 3, 3,
    3, 9, 3,
    3, 3, 9
  ]).reshape((3, 3))

  motif_type = "func"
  num_eigs = 3
  num_clusts = 3

  for sample_weight_type in ["unweighted", "constant", "poisson"]:
    for motif_name in mcut.get_motif_names()[0:15]:
      for mam_weight_type in ["unweighted", "mean", "product"]:
        for type_lap in ["comb", "rw"]:

          # sample a new graph
          adj_mat = mcsa.sample_dsbm(block_sizes, connection_matrix,
                                     weight_matrix, sample_weight_type)

          # run full method
          motif_clust_list = mccl.run_motif_clustering(adj_mat, motif_name,
                                                       motif_type, mam_weight_type,
                                                       "dense", num_eigs,
                                                       type_lap, num_clusts,
                                                       gr_method="dense")

          clusts = motif_clust_list["clusts"]
          comps = motif_clust_list["comps"]

          # answers
          ans_clusts = n * [0] + n * [1] + n * [2]
          ans_clusts = np.take(ans_clusts, comps)

          # score
          ari_score = adjusted_rand_score(clusts, ans_clusts)

          assert ari_score == 1
