#push!(LOAD_PATH,"../src/")

using Documenter
using MotifCluster

makedocs(sitename="MotifCluster.jl",
         modules=MotifCluster,
         pages=["Home" => "index.md",
                "Documentation" => "documentation.md"])

deploydocs(repo="github.com/WGUNDERWOOD/MotifCluster.jl.git")
