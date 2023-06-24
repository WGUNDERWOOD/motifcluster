"""
Get cluster assignments from a spectrum using k-means++.
"""
function cluster_spectrum(vects::Matrix{<:Real}, num_clusts::Int)
    kmeans_plus_plus = kmeans(vects', num_clusts, init = :kmpp,
                              maxiter = 1000, tol = 1e-8)
    cluster_assigns = assignments(kmeans_plus_plus)
    return cluster_assigns
end

"""
Compute the adjusted Rand index between two clusterings.
"""
function adjusted_rand_index(xs::Vector{Int}, ys::Vector{Int})
    @assert length(xs) == length(ys)
    x_vals = sort(unique(xs))
    y_vals = sort(unique(xs))
    n_x_vals = length(x_vals)
    n_y_vals = length(x_vals)
    contingency = zeros(Int, n_x_vals, n_y_vals)
    for i in 1:length(x_vals)
        for j in 1:length(y_vals)
            contingency[i, j] = sum(xs .== x_vals[i] .&& ys .== y_vals[j])
        end
    end
    contingency
    row_sums = sum(contingency, dims=2)
    col_sums = sum(contingency, dims=1)
    a = sum(binomial.(contingency, 2))
    b = sum(binomial.(row_sums, 2)) - a
    c = sum(binomial.(col_sums, 2)) - a
    d = binomial(sum(contingency), 2) - a - b - c
    ari = a - (a + b) * (a + c) / (a + b + c + d)
    ari /= (a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d)
    return ari
end

"""
Run motif-based spectral clustering on the adjacency matrix of a
(weighted directed) network, using a specified motif, motif type,
weighting scheme, embedding dimension, number of clusters and
Laplacian type.
Optionally restrict to the largest connected component before clustering.
"""
function run_motif_clustering(adj_mat::AbstractMatrix{<:Real}, motif_name::String,
                              motif_type::String, mam_weight_type::String, num_eigs::Int,
                              type_lap::String, num_clusts::Int, restrict::Bool)
    @assert motif_type in ["struc", "func"]
    @assert mam_weight_type in ["unweighted", "mean", "product"]
    @assert type_lap in ["comb", "rw"]
    embedding = run_motif_embedding(adj_mat, motif_name, motif_type, mam_weight_type,
                                    num_eigs, type_lap, restrict)
    cluster_assigns = cluster_spectrum(embedding["vects"], num_clusts)
    ans = Dict(
               "adj_mat" => adj_mat,
               "motif_adj_mat" => embedding["motif_adj_mat"],
               "comps" => embedding["comps"],
               "adj_mat_comps" => embedding["adj_mat_comps"],
               "motif_adj_mat_comps" => embedding["motif_adj_mat_comps"],
               "vals" => embedding["vals"],
               "vects" => embedding["vects"],
               "clusts" => cluster_assigns
              )
    return ans
end
