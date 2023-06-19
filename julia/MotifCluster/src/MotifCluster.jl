module MotifCluster

using Clustering
using LinearAlgebra
using SparseArrays
using Graphs

# exports

# includes
include("motifadjacency.jl")
include("utils.jl")
include("indicators.jl")
include("spectral.jl")
include("clustering.jl")
include("sampling.jl")

end
