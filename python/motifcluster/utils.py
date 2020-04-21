import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import diags

def a_b_one(a, b):

  a_sparse = csr_matrix(a)
  b_sparse = csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.T.multiply(b_sparse @ ones_vec)).T

  return csr_matrix(ans)

def a_one_b(a, b):

  a_sparse = csr_matrix(a)
  b_sparse = csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.multiply(ones_vec @ b))

  return csr_matrix(ans)
