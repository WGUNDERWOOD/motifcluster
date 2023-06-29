var documenterSearchIndex = {"docs":
[{"location":"#MotifCluster","page":"Home","title":"MotifCluster","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: motifcluster logo)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A Julia package for motif-based spectral clustering of weighted directed networks.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The MotifCluster package provides implementations of motif-based spectral clustering of weighted directed networks in Julia. These provide the capability for:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Building motif adjacency matrices\nSampling random weighted directed networks\nSpectral embedding with motif adjacency matrices\nMotif-based spectral clustering","category":"page"},{"location":"","page":"Home","title":"Home","text":"The methods are all designed to run quickly on large sparse networks, and are easy to install and use. These methods are based on those described in [Underwood, Elliott and Cucuringu, 2020], which is available at arxiv:2004.01293.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"From the Julia General registry:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add MotifCluster","category":"page"},{"location":"#Dependencies","page":"Home","title":"Dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Aqua\nClustering\nDistributions\nGraphs","category":"page"},{"location":"#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for the MotifCluster package is available on  the web.","category":"page"},{"location":"#Tutorial","page":"Home","title":"Tutorial","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A tutorial for the MotifCluster package is available on Github in the tutorial directory.","category":"page"},{"location":"#Author","page":"Home","title":"Author","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"William George Underwood, Princeton University (maintainer)","category":"page"},{"location":"#Links","page":"Home","title":"Links","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Source code repository on GitHub\nDocumentation on the web.","category":"page"},{"location":"documentation/#Documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"documentation/","page":"Documentation","title":"Documentation","text":"Modules = [MotifCluster]\nPages   = [\"motifadjacency.jl\", \"utils.jl\", \"indicators.jl\", \"spectral.jl\", \"clustering.jl\", \"sampling.jl\"]","category":"page"},{"location":"documentation/#MotifCluster.build_motif_adjacency_matrix-Tuple{AbstractArray{<:Real}, String}","page":"Documentation","title":"MotifCluster.build_motif_adjacency_matrix","text":"Build a motif adjacency matrix from an adjacency matrix. Entry (i, j) of a motif adjacency matrix is the sum of the weights of all motifs containing both nodes i and j.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.a_b_one-Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}","page":"Documentation","title":"MotifCluster.a_b_one","text":"Compute a * (b @ ones) where a, b, ones are square matrices of the same size, and ones contains all entries equal to one. The product * is an entry-wise (Hadamard) product, while @ represents matrix multiplication. This method is more efficient than the naive approach when a or b are sparse.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.a_one_b-Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}","page":"Documentation","title":"MotifCluster.a_one_b","text":"Compute a .* (ones * b) where a, b, ones are square matrices of the same size, and ones contains all entries equal to one.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.dropzeros_killdiag-Tuple{AbstractArray{<:Real}}","page":"Documentation","title":"MotifCluster.dropzeros_killdiag","text":"Set diagonal entries to zero and sparsify.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.get_largest_component-Tuple{AbstractArray{<:Real}}","page":"Documentation","title":"MotifCluster.get_largest_component","text":"Get the indices of the vertices in the largest connected component of a graph from its adjacency matrix.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.get_motif_names-Tuple{}","page":"Documentation","title":"MotifCluster.get_motif_names","text":"Get the names of some common motifs as strings.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_G-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_G","text":"Build the adjacency matrix G.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Gd-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Gd","text":"Build the double-edge adjacency matrix Gd.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Gp-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Gp","text":"Build the product matrix Gp.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Gs-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Gs","text":"Build the single-edge adjacency matrix Gs.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Id-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Id","text":"Build the identity matrix Id.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_J-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_J","text":"Build the directed indicator matrix J.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_J0-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_J0","text":"Build the missing-edge indicator matrix J0.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Jd-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Jd","text":"Build the double-edge indicator matrix Jd.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Je-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Je","text":"Build the edge-and-diagonal matrix Ie.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Jn-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Jn","text":"Build the vertex-distinct indicator matrix Jn.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_Js-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}}","page":"Documentation","title":"MotifCluster.build_Js","text":"Build the single-edge indicator matrix Js.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.build_laplacian-Tuple{AbstractArray{<:Real}, String}","page":"Documentation","title":"MotifCluster.build_laplacian","text":"Build a Laplacian matrix (combinatorial Laplacian or random-walk Laplacian) from a symmetric (weighted) graph adjacency matrix.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.get_first_eigs-Tuple{SparseArrays.SparseMatrixCSC{<:Real, Int64}, Int64}","page":"Documentation","title":"MotifCluster.get_first_eigs","text":"Compute the first few eigenvalues by magnitude and associated eigenvectors of a matrix.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.run_laplace_embedding-Tuple{AbstractArray{<:Real}, Int64, String}","page":"Documentation","title":"MotifCluster.run_laplace_embedding","text":"Run Laplace embedding on a symmetric (weighted) adjacency matrix with a specified number of eigenvalues and eigenvectors.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.run_motif_embedding-Tuple{AbstractArray{<:Real}, String, String, String, Int64, String, Bool}","page":"Documentation","title":"MotifCluster.run_motif_embedding","text":"Calculate a motif adjacency matrix for a given motif and motif type, optionally restrict it to its largest connected component, and then run Laplace embedding with specified Laplacian type and number of eigenvalues and eigenvectors.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.adjusted_rand_index-Tuple{Vector{Int64}, Vector{Int64}}","page":"Documentation","title":"MotifCluster.adjusted_rand_index","text":"Compute the adjusted Rand index between two clusterings.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.cluster_spectrum-Tuple{Matrix{<:Real}, Int64}","page":"Documentation","title":"MotifCluster.cluster_spectrum","text":"Get cluster assignments from a spectrum using k-means++.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.run_motif_clustering-Tuple{AbstractMatrix{<:Real}, String}","page":"Documentation","title":"MotifCluster.run_motif_clustering","text":"Run motif-based spectral clustering on the adjacency matrix of a (weighted directed) network, using a specified motif, motif type, weighting scheme, embedding dimension, number of clusters and Laplacian type. Optionally restrict to the largest connected component before clustering.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.demonstration_graph-Tuple{}","page":"Documentation","title":"MotifCluster.demonstration_graph","text":"Generate a small graph for demonstrations.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.random_sparse_matrix-Tuple{Int64, Int64, Float64, String, Real}","page":"Documentation","title":"MotifCluster.random_sparse_matrix","text":"Build a sparse matrix of size m, n with non-zero probability p. Edge weights can be unweighted, constant-weighted or Poisson-weighted.\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.sample_bsbm-Tuple{Vector{Int64}, Vector{Int64}, Matrix{<:Real}}","page":"Documentation","title":"MotifCluster.sample_bsbm","text":"Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic block model (BSBM).\n\n\n\n\n\n","category":"method"},{"location":"documentation/#MotifCluster.sample_dsbm-Tuple{Vector{Int64}, Matrix{<:Real}}","page":"Documentation","title":"MotifCluster.sample_dsbm","text":"Sample the (weighted) adjacency matrix of a (weighted) directed stochastic block model (DSBM) with specified parameters.\n\n\n\n\n\n","category":"method"}]
}
