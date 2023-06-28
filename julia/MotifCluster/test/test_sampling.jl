@testset verbose = true "Clustering" begin
    @testset verbose = true "sample_dsbm_unweighted" begin
        Random.seed!(9349)
        sample_weight_type = "unweighted"
        block_sizes = [2, 3]
        connection_matrix = [0.4 0.5; 0.6 0.7]
        n_reps = 2000
        weight_matrix = nothing
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                     weight_matrix=weight_matrix,
                                     sample_weight_type=sample_weight_type)
        G_vals = [0, 1]
        @test all(g in G_vals for g in G)
        G = zeros(5, 5)
        for rep in 1:n_reps
            G += MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                          weight_matrix=weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 0.4 0.5 0.5 0.5;
                      0.4 0 0.5 0.5 0.5;
                      0.6 0.6 0 0.7 0.7;
                      0.6 0.6 0.7 0 0.7;
                      0.6 0.6 0.7 0.7 0])
        @test isapprox(G, ans, atol=0.05)
    end

    @testset verbose = true "sample_dsbm_constant_weighted" begin
        Random.seed!(2839)
        sample_weight_type = "constant"
        block_sizes = [2, 3]
        connection_matrix = [0.4 0.5; 0.6 0.7]
        n_reps = 2000
        weight_matrix = [20 30; 40 50]
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                     weight_matrix=weight_matrix,
                                     sample_weight_type=sample_weight_type)
        G_vals = [0, 20, 30, 40, 50]
        @test all(g in G_vals for g in G)

        G = zeros(5, 5)
        for rep in 1:n_reps
            G += MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                          weight_matrix=weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 8 15 15 15;
                      8 0 15 15 15;
                      24 24 0 35 35;
                      24 24 35 0 35;
                      24 24 35 35 0])
        @test isapprox(G, ans, atol=3)
    end

    @testset verbose = true "sample_dsbm_poisson_weighted" begin
        Random.seed!(2838)
        sample_weight_type = "poisson"
        block_sizes = [2, 3]
        connection_matrix = [0.4 0.5; 0.6 0.7]
        n_reps = 2000
        weight_matrix = [20 30; 40 50]
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                     weight_matrix=weight_matrix,
                                     sample_weight_type=sample_weight_type)
        @test all(G .== floor.(G))
        @test all(G .>= 0)
        G = zeros(5, 5)
        for rep in 1:n_reps
            G += MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                          weight_matrix=weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 8 15 15 15;
                      8 0 15 15 15;
                      24 24 0 35 35;
                      24 24 35 0 35;
                      24 24 35 35 0])
        @test isapprox(G, ans, atol=3)
    end

    @testset verbose = true "sample_dsbm_large" begin
        Random.seed!(2238)
        n = Int(1e5)
        block_sizes = [n]
        connection_matrix = reshape([10 / n], (1, 1))
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix;
                                     weight_matrix=nothing,
                                     sample_weight_type="unweighted")
        @test size(G) == (n, n)
    end

    @testset verbose = true "sample_bsbm_unweighted" begin
        Random.seed!(9424)
        sample_weight_type = "unweighted"
        source_block_sizes = [1, 2]
        dest_block_sizes = [1, 1, 1]
        bipartite_connection_matrix = [0.3 0.4 0.5; 0.6 0.7 0.8]
        n_reps = 2000
        bipartite_weight_matrix = nothing
        G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                     bipartite_connection_matrix;
                                     bipartite_weight_matrix=bipartite_weight_matrix,
                                     sample_weight_type=sample_weight_type)
        G_vals = [0, 1]
        @test all(g in G_vals for g in G)
        G = zeros(6, 6)
        for rep in 1:n_reps
            G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                          bipartite_connection_matrix;
                                          bipartite_weight_matrix=bipartite_weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 0 0 0.3 0.4 0.5;
                      0 0 0 0.6 0.7 0.8;
                      0 0 0 0.6 0.7 0.8;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0])
        @test isapprox(G, ans, atol=0.05)
    end

    @testset verbose = true "sample_bsbm_constant_weighted" begin
        Random.seed!(7482)
        sample_weight_type = "constant"
        source_block_sizes = [1, 2]
        dest_block_sizes = [1, 1, 1]
        bipartite_connection_matrix = [0.3 0.4 0.5; 0.6 0.7 0.8]
        n_reps = 2000
        bipartite_weight_matrix = [10 20 30; 40 50 60]
        G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                     bipartite_connection_matrix;
                                     bipartite_weight_matrix=bipartite_weight_matrix,
                                     sample_weight_type=sample_weight_type)
        G_vals = [0, 10, 20, 30, 40, 50, 60]
        @test all(g in G_vals for g in G)
        G = zeros(6, 6)
        for rep in 1:n_reps
            G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                          bipartite_connection_matrix;
                                          bipartite_weight_matrix=bipartite_weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 0 0 3 8 15;
                      0 0 0 24 35 48;
                      0 0 0 24 35 48;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0])
        @test isapprox(G, ans, atol=3)
    end

    @testset verbose = true "sample_bsbm_poisson_weighted" begin
        Random.seed!(7482)
        sample_weight_type = "poisson"
        source_block_sizes = [1, 2]
        dest_block_sizes = [1, 1, 1]
        bipartite_connection_matrix = [0.3 0.4 0.5; 0.6 0.7 0.8]
        n_reps = 2000
        bipartite_weight_matrix = [10 20 30; 40 50 60]
        G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                     bipartite_connection_matrix;
                                     bipartite_weight_matrix=bipartite_weight_matrix,
                                     sample_weight_type=sample_weight_type)
        @test all(G .== floor.(G))
        @test all(G .>= 0)
        G = zeros(6, 6)
        for rep in 1:n_reps
            G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes,
                                          bipartite_connection_matrix;
                                          bipartite_weight_matrix=bipartite_weight_matrix,
                                          sample_weight_type=sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0 0 0 3 8 15;
                      0 0 0 24 35 48;
                      0 0 0 24 35 48;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0;
                      0 0 0 0 0 0])
        @test isapprox(G, ans, atol=3)
    end
end
