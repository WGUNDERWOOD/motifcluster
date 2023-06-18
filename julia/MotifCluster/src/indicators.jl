"""
Functions for building adjacency and indicator matrices
are in `motifcluster.indicators`.
"""

"""
Build the adjacency matrix `G`.
"""
function build_G(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = dropzeros_killdiag(adj_mat)
    return G
end

"""
Build the directed indicator matrix `J`.
"""
function build_J(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = build_G(adj_mat)
    J = dropzeros_killdiag(1 .* (G .> 0))
    return J
end

"""
Build the single-edge adjacency matrix `Gs`.
"""
function build_Gs(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = build_G(adj_mat)
    J = build_J(adj_mat)
    Gs = dropzeros_killdiag(G - G .* J')
    return Gs
end

"""
Build the single-edge indicator matrix `Js`.
"""
function build_Js(adj_mat::SparseMatrixCSC{<:Real, Int})
    Gs = build_Gs(adj_mat)
    Js = dropzeros_killdiag(1 .* (Gs .> 0))
    return Js
end

"""
Build the double-edge adjacency matrix `Gd`.
"""
function build_Gd(adj_mat::SparseMatrixCSC{<:Real, Int})
    J = build_J(adj_mat)
    G = build_G(adj_mat)
    Gd = dropzeros_killdiag((G .+ G') .* J .* J')
    return Gd
end

"""
Build the double-edge indicator matrix `Jd`.
"""
function build_Jd(adj_mat::SparseMatrixCSC{<:Real, Int})
    Gd = build_Gd(adj_mat)
    Jd = dropzeros_killdiag(1 .* (Gd .> 0))
    return Jd
end

"""
Build the missing-edge indicator matrix `J0`.
"""
function build_J0(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = build_G(adj_mat)
    J0 = dropzeros_killdiag((G .+ G') .== 0)
    return J0
end

"""
Build the vertex-distinct indicator matrix `Jn`.
"""
function build_Jn(adj_mat::SparseMatrixCSC{<:Real, Int})
    Jn = dropzeros_killdiag(ones(size(adj_mat)))
    return Jn
end

"""
Build the identity matrix `Id`.
"""
function build_Id(adj_mat::SparseMatrixCSC{<:Real, Int})
    Id = sparse(Matrix(I, size(adj_mat)))
    return Id
end

"""
Build the edge-and-diagonal matrix `Ie`.
"""
function build_Je(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = build_G(adj_mat)
    Id = build_Id(adj_mat)
    Je = Id .+ ((G .+ G') .> 0)
    return Je
end

"""
Build the product matrix `Gp`.
"""
function build_Gp(adj_mat::SparseMatrixCSC{<:Real, Int})
    G = build_G(adj_mat)
    Gp = dropzeros_killdiag(G .* G')
    return Gp
end
