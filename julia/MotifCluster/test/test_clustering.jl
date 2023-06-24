@testset verbose = true "Clustering" begin

    @testset verbose = true "cluster_spectrum" begin

        Random.seed!(1)
        n = 10
        vects_1 = rand(Normal(), (n, 3))
        vects_2 = hcat(rand(Normal(), (n, 1)), rand(Normal(4), (n, 2)))

        vals = Int[]
        vects =  vcat(vects_1, vects_2)

        clust_ans = [[1 for _ in 1:n]; [2 for _ in 1:n]]
        clust = MotifCluster.cluster_spectrum(vects, 2)

        if clust[1] == 2
            clust .= 3 .- clust
        end

        @test clust == clust_ans
    end

    @testset verbose = true "run_motif_clustering" begin
        Random.seed!(3965)
        n = 50
        block_sizes = [n, n, n]
        connection_matrix = [0.9 0.2 0.2; 0.2 0.9 0.2; 0.2 0.2 0.9]
        weight_matrix = [9 2 2; 2 9 2; 2 2 9]
        motif_type = "func"
        num_eigs = 3
        num_clusts = 3
        for sample_weight_type in ["unweighted", "constant", "poisson"]
            for motif_name in MotifCluster.get_motif_names()[1:15]
                for mam_weight_type in ["unweighted", "mean", "product"]
                    for type_lap in ["comb", "rw"]
                        adj_mat_sparse = MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                                                  weight_matrix = weight_matrix,
                                                                  sample_weight_type = sample_weight_type)
                        adj_mat_dense = Matrix(adj_mat_sparse)
                        for adj_mat in [adj_mat_sparse, adj_mat_dense]
                            motif_clust_list = MotifCluster.run_motif_clustering(adj_mat, motif_name;
                                                                                 motif_type = motif_type,
                                                                                 mam_weight_type = mam_weight_type,
                                                                                 num_eigs = num_eigs,
                                                                                 type_lap = type_lap,
                                                                                 num_clusts = num_clusts,
                                                                                 restrict = true)
                            clusts = motif_clust_list["clusts"]
                            comps = motif_clust_list["comps"]
                            ans_clusts = [repeat([1], n); repeat([2], n); repeat([3], n)]
                            ans_clusts = ans_clusts[comps]
                            ari_score = MotifCluster.adjusted_rand_index(clusts, ans_clusts)
                            if ari_score < 1
                                println(motif_name)
                                println(sample_weight_type)
                                println(type_lap)
                            end
                            @test ari_score == 1
                        end
                    end
                end
            end
        end
    end
end
