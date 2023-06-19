"""
Build a sparse matrix of size `m, n` with non-zero probability `p`.
Edge weights can be unweighted, constant-weighted or
Poisson-weighted.
"""
function random_sparse_matrix(m::Int, n::Int, p::Float64,
        sample_weight_type::String, w::Int)
    if sample_weight_type == "constant"
        return w * (sprand(m, n, p) .> 0)
    elseif sample_weight_type == "poisson"
        return sprandn(Distributions.Poisson(w), m, n, p)
    else
        return sprand(Bool, m, n, p)
    end
end

"""
Sample the (weighted) adjacency matrix of a (weighted) directed stochastic
block model (DSBM) with specified parameters.
"""
function sample_dsbm(block_sizes::Vector{Int}, connection_matrix::Matrix{Float64},
        weight_matrix::Union{Matrix{Float64}, Nothing}, sample_weight_type::String)

    # check args
    @assert all(block_sizes .> 0)
    @assert length(block_sizes) == size(connection_matrix, 1) == size(connection_matrix, 2)
    @assert all(0 .<= connection_matrix .<= 1)

    if sample_weight_type != "unweighted"
        @assert weight_matrix != nothing
        @assert length(block_sizes) == size(weight_matrix, 1) == size(weight_matrix, 2)
        @assert all(weight_matrix .>= 0)
    end

    # initialize variables
    k = length(block_sizes)
    block_list = SparseMatrixCSC{<:Real, Int}[]

    for i in 1:k
        for j in 1:k

            # block parameters
            ni = block_sizes[i]
            nj = block_sizes[j]
            p = connection_matrix[i, j]

            # generate block
            if sample_weight_type == "unweighted"
                w = 1
            else
                w = weight_matrix[i, j]
            end

            block = random_sparse_matrix(ni, nj, p, sample_weight_type, w)
            push!(block_list, block)
        end
    end
    adj_mat = hvcat(size(connection_matrix, 2), block_list...)
    adj_mat = dropzeros_killdiag(adj_mat)
    return adj_mat
end

#=
function sample_bsbm(source_block_sizes, dest_block_sizes,
                bipartite_connection_matrix,
                bipartite_weight_matrix=None,
                sample_weight_type="unweighted")

  """
  Sample a bipartite stochastic block model (BSBM).

  Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic
  block model (BSBM) with specified parameters.

  Parameters
  ----------
  source_block_sizes : list of int
    A list containing the size of each block of source vertices.
  dest_block_sizes : list of int
    A list containing the size of each block of destination vertices.
  bipartite_connection_matrix : matrix
    A matrix containing the source block to destination block
    connection probabilities.
  bipartite_weight_matrix : matrix
    A matrix containing the source block to destination block weight
    parameters. Unused for `sample_weight_type = "constant"`.
    Defaults to `None`.
  sample_weight_type : str
    The type of weighting scheme.
    One of `"unweighted"`, `"constant"` or `"poisson"`.

  Returns
  -------
  adj_mat : sparse matrix
    A randomly sampled (weighted) adjacency matrix of a BSBM.

  Examples
  --------
  >>> source_block_sizes = [10, 10]
  >>> dest_block_sizes = [10, 10, 10]
  >>> bipartite_connection_matrix = np.array([0.8, 0.5, 0.1, 0.1, 0.5, 0.8]).reshape((2, 3))
  >>> bipartite_weight_matrix = np.array([20, 10, 2, 2, 10, 20]).reshape((2, 3))
  >>> sample_bsbm(block_sizes, bipartite_connection_matrix,
  ...   bipartite_weight_matrix, "poisson")
  """

  # check args
  assert source_block_sizes == [int(x) for x in source_block_sizes]
  assert dest_block_sizes == [int(x) for x in dest_block_sizes]
  assert all(x > 0 for x in source_block_sizes)
  assert all(x > 0 for x in dest_block_sizes)
  assert len(source_block_sizes) == bipartite_connection_matrix.shape[0]
  assert len(dest_block_sizes) == bipartite_connection_matrix.shape[1]
  assert (bipartite_connection_matrix >= 0).all()
  assert (bipartite_connection_matrix <= 1).all()
  assert sample_weight_type in ["unweighted", "constant", "poisson"]

  if sample_weight_type != "unweighted":
    assert bipartite_weight_matrix is not None
    assert len(source_block_sizes) == bipartite_weight_matrix.shape[0]
    assert len(dest_block_sizes) == bipartite_weight_matrix.shape[1]
    assert (bipartite_weight_matrix >= 0).all()

  # initialize parameters
  ks = len(source_block_sizes)
  kd = len(dest_block_sizes)
  zeros_ss = np.zeros((ks, ks))
  zeros_d = np.zeros((kd, ks + kd))

  # build block sizes vector
  block_sizes = source_block_sizes + dest_block_sizes

  # build connection matrix
  connection_matrix = np.block([[zeros_ss, bipartite_connection_matrix], [zeros_d]])

  # build weight matrix
  if bipartite_weight_matrix is not None:
    weight_matrix = np.block([[zeros_ss, bipartite_weight_matrix], [zeros_d]])

  else:
    weight_matrix = None

  # sample BSBM
  adj_mat = sample_dsbm(block_sizes, connection_matrix,
                        weight_matrix, sample_weight_type)

  return adj_mat

  =#

"""
Generate a small graph for demonstrations.
"""
function demonstration_graph()

    adj_mat = sparse([0 0  0  0  0  0  0  0  0  0  0 0;
                      2 0  3  0  6  8  0  0  0  0  0 0;
                      0 0  0  0  0 10  0  0  0  0  0 0;
                      0 0  0  0  0  0  0  0  0  0  0 0;
                      4 5  0  0  0 14  0  0 18 19  0 0;
                      0 7  9  0 13  0  0  0  0  0 21 0;
                      0 0 11 12  0 15  0 17  0  0  0 0;
                      0 0  0  0  0  0 16  0  0  0  0 0;
                      0 0  0  0  0  0  0  0  0 24  0 0;
                      0 0  0  0  0 20  0  0  0  0  0 0;
                      0 0  0  0  0  0 22  0  0  0  0 0;
                      0 0  0  0  0  0 23  0  0  0  0 0])

    return adj_mat
  end
