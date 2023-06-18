#=
def _a_b_one(a_mat, b_mat):

  """
  Compute a right-multiplication with the ones matrix.

  Compute `a * (b @ one_mat)` where `a`, `b`,
  `ones_mat` are square matrices of the same size,
  and `ones_mat` contains all entries equal to one.
  The product `*` is an entry-wise (Hadamard) product,
  while `@` represents matrix multiplication.
  This method is more efficient than the naive approach
  when `a` or `b` are sparse.

  Parameters
  ----------
  a, b : matrix
    Square matrices of the same size.

  Returns
  -------
  sparse matrix
    The sparse square matrix `a * (b @ one_mat)`.
  """

  if not sparse.issparse(a_mat):
    a_mat = sparse.csr_matrix(a_mat)

  if not sparse.issparse(b_mat):
    b_mat = sparse.csr_matrix(b_mat)

  n = a_mat.shape[0]
  ones_vec = np.ones(n)
  ans = (a_mat.T.multiply(b_mat @ ones_vec)).T

  return ans


def _a_one_b(a_mat, b_mat):

  """
  Compute a left-multiplication with the ones matrix.

  Compute `a * (one_mat @ b)` where `a`, `b`,
  `ones_mat` are square matrices of the same size,
  and `ones_mat` contains all entries equal to one.
  The product `*` is an entry-wise (Hadamard) product,
  while `@` represents matrix multiplication.
  This method is more efficient than the naive approach
  when `a` or `b` are sparse.

  Parameters
  ----------
  a, b : matrix
    Square matrices of the same size.

  Returns
  -------
  sparse matrix
    The sparse square matrix `a * (one_mat @ b)`.
  """

  if not sparse.issparse(a_mat):
    a_mat = sparse.csr_matrix(a_mat)

  if not sparse.issparse(b_mat):
    b_mat = sparse.csr_matrix(b_mat)

  n = a_mat.shape[0]
  ones_vec = np.ones(n)
  ans = (a_mat.multiply(ones_vec @ b_mat))

  return ans
  =#


"""
Set diagonal entries to zero and sparsify.
"""
function dropzeros_killdiag(some_mat::AbstractArray{<:Real})::SparseMatrixCSC{Float64, Int}
    I = Diagonal(ones(size(some_mat, 1)))
    ans = dropzeros(sparse(some_mat .- some_mat .* I))
    return ans
end

  #=

def get_largest_component(adj_mat, gr_method):

  """
  Get largest connected component.

  Get the indices of the vertices in the largest connected
  component of a graph from its adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    An adjacency matrix of a graph.
  gr_method : str
    Format to use before building the graph.
    One of `"sparse"` or `"dense"`.

  Returns
  -------
  verts_to_keep : list
    A list of indices corresponding to the vertices in the largest
    connected component.

  Examples
  --------
  >>> adj_mat = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0]).reshape((3, 3))
  >>> get_largest_component(adj_mat)
  """

  if gr_method == "sparse":

    if not sparse.issparse(adj_mat):
      adj_mat = sparse.csr_matrix(adj_mat)

    gr = nx.from_scipy_sparse_array(adj_mat > 0)

  else:

    if isinstance(adj_mat, np.ndarray): # pylint: disable=else-if-used
      gr = nx.from_numpy_array(1 * np.array(adj_mat > 0))

    else:
      gr = nx.from_numpy_array(1 * (adj_mat > 0).toarray())

  verts_to_keep = max(nx.connected_components(gr), key=len)
  verts_to_keep = sorted(verts_to_keep)

  return verts_to_keep

  =#

"""
Get the names of some common motifs as strings.
"""
function get_motif_names()
    motif_names = ["Ms", "Md"]
    for i in 1:14
        motif_name = "M" * string(i)
        push!(motif_names, motif_name)
    end
    push!(motif_names, "Mcoll")
    push!(motif_names, "Mexpa")
    return motif_names
end
