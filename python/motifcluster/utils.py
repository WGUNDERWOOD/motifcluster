import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import diags

def a_b_one(a, b):

  n = np.shape(a)[0]
  ones_vec = np.ones((n, 1))
  diagonal = (b @ ones_vec).reshape((n))
  ans = diags(diagonal) @ a

  return csr_matrix(ans)
