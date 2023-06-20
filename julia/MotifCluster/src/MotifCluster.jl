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

# includes
include("motifadjacency.jl")
include("utils.jl")
include("indicators.jl")
include("spectral.jl")
include("clustering.jl")
include("sampling.jl")

end
