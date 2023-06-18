struct Spectrum
    vals::Vector{Float64}
    vects::Matrix{Float64}
end

"""
Compute the first few eigenvalues by magnitude and
associated eigenvectors of a matrix.
"""
function get_first_eigs(some_mat::SparseMatrixCSC{<:Real, Int}, num_eigs::Int)
    @assert num_eigs >= 1
    spectrum = eigen(Matrix(some_mat))
    vals = spectrum.values[1:num_eigs]
    vects = spectrum.vectors[:, 1:num_eigs]
    return Spectrum(vals, vects)
end

"""
Build a Laplacian matrix (combinatorial Laplacian or random-walk Laplacian)
from a symmetric (weighted) graph adjacency matrix.
"""
function build_laplacian(adj_mat::SparseMatrixCSC{<:Real, Int}, type_lap::String)
    degs = vec(sum(adj_mat, dims=1))
    if type_lap == "comb"
        degs_mat = sparse(diagm(degs))
        L = degs_mat - adj_mat
    elseif type_lap == "rw"
        @assert all(degs .> 0)
        inv_degs_mat = sparse(diagm(1 ./ degs))
        L = sparse(I - inv_degs_mat * adj_mat)
    end
    return L
end

"""
Run Laplace embedding on a symmetric (weighted) adjacency matrix
with a specified number of eigenvalues and eigenvectors.
"""
function run_laplace_embedding(adj_mat::SparseMatrixCSC, num_eigs::Int, type_lap::String)
    @assert num_eigs >= 1
    laplacian = build_laplacian(adj_mat, type_lap)
    spectrum = get_first_eigs(laplacian, num_eigs)
    return spectrum
end

"""
Calculate a motif adjacency matrix for a given motif and motif type,
optionally restrict it to its largest connected component,
and then run Laplace embedding with specified Laplacian type and
number of eigenvalues and eigenvectors.
"""
function run_motif_embedding(adj_mat::AbstractArray{<:Real}, motif_name::String,
        motif_type::String, mam_weight_type::String,
        num_eigs::Int, type_lap::String, restrict::Bool)
    @assert num_eigs >= 1
    motif_adj_mat = build_motif_adjacency_matrix(adj_mat, motif_name, motif_type, mam_weight_type)
    if restrict
        comps = get_largest_component(motif_adj_mat)
        adj_mat_comps = sparse(adj_mat[comps, comps])
        motif_adj_mat_comps = motif_adj_mat[comps, comps]
        spect = run_laplace_embedding(motif_adj_mat_comps, num_eigs, type_lap)
    else
        comps = nothing
        adj_mat_comps = nothing
        motif_adj_mat_comps = nothing
        spect = run_laplace_embedding(motif_adj_mat, num_eigs, type_lap)
    end
    embedding = MotifEmbedding(adj_mat, motif_adj_mat, comps, adj_mat_comps,
                               motif_adj_mat_comps, spect.vals, spect.vects)
    return embedding
end
