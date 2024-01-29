using MotifCluster
using Test
using Random
using Distributions
using SparseArrays
using LinearAlgebra
using Aqua

@testset verbose = true "MotifCluster" begin
    Aqua.test_ambiguities(MotifCluster)
    Aqua.test_unbound_args(MotifCluster)
    Aqua.test_undefined_exports(MotifCluster)
    Aqua.test_project_extras(MotifCluster)
    Aqua.test_stale_deps(MotifCluster, ignore=[:Aqua, :Documenter])
    Aqua.test_deps_compat(MotifCluster)
    include("test_clustering.jl")
    include("test_indicators.jl")
    include("test_spectral.jl")
    include("test_motifadjacency.jl")
    include("test_utils.jl")
    include("test_sampling.jl")
end
