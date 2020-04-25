from motifcluster import clustering as mccl
import numpy as np
from numpy import random as rd

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

  print(clust_ans)

  assert (clust == clust_ans).all()
