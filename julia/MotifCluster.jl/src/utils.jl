"""
Compute `a * (b @ ones)` where `a`, `b`,
`ones` are square matrices of the same size,
and `ones` contains all entries equal to one.
The product `*` is an entry-wise (Hadamard) product,
while `@` represents matrix multiplication.
This method is more efficient than the naive approach
when `a` or `b` are sparse.
"""
function a_b_one(a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real})
    ones_vec = ones(size(a, 1), 1)
    ans = (a .* (b * ones_vec))
    return ans
end

"""
Compute `a .* (ones * b)` where `a`, `b`,
`ones` are square matrices of the same size,
and `ones` contains all entries equal to one.
"""
function a_one_b(a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real})
    ones_vec = ones(size(a, 1), 1)
    ans = (a .* (ones_vec' * b))
    return ans
end

"""
Set diagonal entries to zero and sparsify.
"""
function dropzeros_killdiag(some_mat::AbstractArray{<:Real})::SparseMatrixCSC{Float64,Int}
    I = Diagonal(ones(size(some_mat, 1)))
    ans = dropzeros(sparse(some_mat .- some_mat .* I))
    return ans
end

"""
Get the indices of the vertices in the largest connected
component of a graph from its adjacency matrix.
"""
function get_largest_component(adj_mat::AbstractArray{<:Real})
    gr = Graph(adj_mat + adj_mat' .> 0)
    comps = Graphs.connected_components(gr)
    max_size = maximum(length(c) for c in comps)
    verts_to_keep = [c for c in comps if length(c) == max_size][]
    return sort(verts_to_keep)
end

"""
Get the names of some common motifs as strings.
"""
function get_motif_names()
    motif_names = ["Ms", "Md"]
    for i in 1:13
        motif_name = "M" * string(i)
        push!(motif_names, motif_name)
    end
    push!(motif_names, "Mcoll")
    push!(motif_names, "Mexpa")
    return motif_names
end
