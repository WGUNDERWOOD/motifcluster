module MotifCluster

using Clustering
using LinearAlgebra
using SparseArrays
using Graphs
using Distributions

# exports
export build_motif_adjacency_matrix
export sample_dsbm
export sample_bsbm
export get_motif_names
export run_motif_clustering
export adjusted_rand_index

# includes
include("motifadjacency.jl")
include("utils.jl")
include("indicators.jl")
include("spectral.jl")
include("clustering.jl")
include("sampling.jl")

end
