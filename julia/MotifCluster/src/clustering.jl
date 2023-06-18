"""
Get cluster assignments from a spectrum using k-means++.
"""
function cluster_spectrum(spectrum::Spectrum, num_clusts::Int)
    vects = spectrum.vects[:, 1:end]
    kmeans_plus_plus = kmeans(vects', num_clusts, init = :kmpp)
    cluster_assigns = assignments(kmeans_plus_plus)
    return cluster_assigns
end


"""
Run motif-based spectral clustering on the adjacency matrix of a
(weighted directed) network, using a specified motif, motif type,
weighting scheme, embedding dimension, number of clusters and
Laplacian type.
Optionally restrict to the largest connected component before clustering.
"""
function run_motif_clustering(adj_mat::Matrix{Float64}, motif_name::String,
                              motif_type::String="struc",
                              mam_weight_type::String="unweighted",
                              mam_method::String="sparse",
                              num_eigs::Int=2,
                              type_lap::String="comb",
                              num_clusts::Int=2,
                              restrict::Bool=true,
                              gr_method::String="sparse")

    @assert motif_type in ["struc", "func"]
    @assert mam_weight_type in ["unweighted", "mean", "product"]
    @assert mam_method in ["sparse", "dense"]
    @assert type_lap in ["comb", "rw"]

    spectrum = run_motif_embedding(adj_mat, motif_name, motif_type, mam_weight_type,
                                   mam_method, num_eigs, type_lap, restrict, gr_method)

    cluster_assigns = cluster_spectrum(spectrum, num_clusts)

    ans = Dict(
               "adj_mat" => adj_mat,
               "motif_adj_mat" => spectrum["motif_adj_mat"],
               "comps" => spectrum["comps"],
               "adj_mat_comps" => spectrum["adj_mat_comps"],
               "motif_adj_mat_comps" => spectrum["motif_adj_mat_comps"],
               "vals" => spectrum["vals"],
               "vects" => spectrum["vects"],
               "clusts" => cluster_assigns
              )

    return ans
end
