import numpy as np
import networkx as nx
from scipy import sparse

def a_b_one(a, b):

  a_sparse = sparse.csr_matrix(a)
  b_sparse = sparse.csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.T.multiply(b_sparse @ ones_vec)).T

  return sparse.csr_matrix(ans)


def a_one_b(a, b):

  a_sparse = sparse.csr_matrix(a)
  b_sparse = sparse.csr_matrix(b)
  n = a.shape[0]
  ones_vec = np.ones(n)
  ans = (a_sparse.multiply(ones_vec @ b))

  return sparse.csr_matrix(ans)


def drop0_killdiag(some_mat):

  ans = sparse.csr_matrix(some_mat)
  I = sparse.identity(ans.shape[0])
  ans = ans - I.multiply(ans)

  return(sparse.csr_matrix(ans))


def get_largest_component(adj_mat):

  adj_mat_sparse = sparse.csr_matrix(adj_mat)
  gr = nx.from_scipy_sparse_matrix(adj_mat_sparse > 0)
  verts_to_keep = max(nx.connected_components(gr), key=len)
  verts_to_keep = sorted(verts_to_keep)

  return verts_to_keep
