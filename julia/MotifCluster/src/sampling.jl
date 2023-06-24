"""
Build a sparse matrix of size `m, n` with non-zero probability `p`.
Edge weights can be unweighted, constant-weighted or
Poisson-weighted.
"""
function random_sparse_matrix(m::Int, n::Int, p::Float64, sample_weight_type::String, w::Real)
    if sample_weight_type == "constant"
        return w * sprand(Bool, m, n, p)
    elseif sample_weight_type == "poisson"
        d = Distributions.Poisson(w)
        return quantile.(d, sprand(m, n, p))
    else
        return sprand(Bool, m, n, p)
    end
end

"""
Sample the (weighted) adjacency matrix of a (weighted) directed stochastic
block model (DSBM) with specified parameters.
"""
function sample_dsbm(block_sizes::Vector{Int}, connection_matrix::Matrix{<:Real};
        weight_matrix::Union{Matrix{<:Real}, Nothing} = nothing,
        sample_weight_type::String = "unweighted")

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

"""
Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic block model (BSBM).
"""
function sample_bsbm(source_block_sizes::Vector{Int}, dest_block_sizes::Vector{Int},
        bipartite_connection_matrix::Matrix{<:Real};
        bipartite_weight_matrix::Union{Matrix{<:Real}, Nothing} = nothing,
        sample_weight_type::String = "unweighted")

    # check args
    @assert all(source_block_sizes .> 0)
    @assert all(dest_block_sizes .> 0)
    @assert length(source_block_sizes) == size(bipartite_connection_matrix, 1)
    @assert length(dest_block_sizes) == size(bipartite_connection_matrix, 2)
    @assert all(0 .<= bipartite_connection_matrix .<= 1)

    if sample_weight_type != "unweighted"
        @assert !isnothing(bipartite_weight_matrix)
        @assert length(source_block_sizes) == size(bipartite_weight_matrix, 1)
        @assert length(dest_block_sizes) == size(bipartite_weight_matrix, 2)
        @assert all(bipartite_weight_matrix .>= 0)
    end

    # initialize parameters
    ks = length(source_block_sizes)
    kd = length(dest_block_sizes)
    zeros_ss = zeros(ks, ks)
    zeros_d = zeros(kd, ks + kd)

    # build block sizes vector
    block_sizes = [source_block_sizes; dest_block_sizes]

    # build connection matrix
    connection_matrix = hvcat((2, 1), zeros_ss, bipartite_connection_matrix, zeros_d)

    # build weight matrix
    if !isnothing(bipartite_weight_matrix)
        weight_matrix = hvcat((2, 1), zeros_ss, bipartite_weight_matrix, zeros_d)
    else
        weight_matrix = nothing
    end

    # sample BSBM
    adj_mat = sample_dsbm(block_sizes, connection_matrix;
                          weight_matrix = weight_matrix,
                          sample_weight_type = sample_weight_type)
    return adj_mat
end


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
