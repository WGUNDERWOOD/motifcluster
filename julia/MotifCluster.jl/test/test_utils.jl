@testset verbose = true "Utils" begin
    @testset verbose = true "a_b_one" begin
        a = sparse(reshape(collect(-4:4), (3, 3)))'
        b = sparse(reshape(collect(-1:7), (3, 3)))'
        ans = sparse([0 0 0; -9 0 9; 36 54 72])
        @test isapprox(MotifCluster.a_b_one(a, b), ans)
    end

    @testset verbose = true "a_one_b" begin
        a = sparse(reshape(collect(-4:4), (3, 3)))'
        b = sparse(reshape(collect(-1:7), (3, 3)))'
        ans = sparse([-24 -27 -24; -6 0 12; 12 27 48])
        @test isapprox(MotifCluster.a_one_b(a, b), ans)
    end

    @testset verbose = true "dropzeros_killdiag" begin
        adj_mat = sparse(reshape(collect(-1:7), (3, 3)))'
        ans = sparse([0 0 1; 2 0 4; 5 6 0])
        @test isapprox(MotifCluster.dropzeros_killdiag(adj_mat), ans)
    end

    @testset verbose = true "get_largest_component" begin
        adj_mat_dense = [0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 2; 3 0 0 0 0; 0 4 0 0 0]
        adj_mat_sparse = sparse(adj_mat_dense)
        ans = [2, 3, 5]
        @test MotifCluster.get_largest_component(adj_mat_sparse) == ans
        @test MotifCluster.get_largest_component(adj_mat_dense) == ans
    end
end
