#push!(LOAD_PATH,"../src/")

using Documenter
using MotifCluster

makedocs(sitename="MotifCluster.jl",
         modules=MotifCluster,
         remotes=nothing,
         repo = Remotes.GitHub("WGUNDERWOOD", "motifcluster"),
         format = Documenter.HTML(edit_link = "main",
                                  repolink = "https://github.com/WGUNDERWOOD/motifcluster"),
         pages=["Home" => "index.md",
                "Documentation" => "documentation.md"])

deploydocs(repo="github.com/WGUNDERWOOD/motifcluster.git")
