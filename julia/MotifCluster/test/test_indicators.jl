@testset verbose = true "Indicators" begin

    adj_mat = [0 2 0; 3 0 4; 0 0 0]

    G  = sparse([0 2 0; 3 0 4; 0 0 0])
    J  = sparse([0 1 0; 1 0 1; 0 0 0])
    Gs = sparse([0 0 0; 0 0 4; 0 0 0])
    Js = sparse([0 0 0; 0 0 1; 0 0 0])
    Gd = sparse([0 5 0; 5 0 0; 0 0 0])
    Jd = sparse([0 1 0; 1 0 0; 0 0 0])
    J0 = sparse([0 0 1; 0 0 0; 1 0 0])
    Jn = sparse([0 1 1; 1 0 1; 1 1 0])
    Id = sparse([1 0 0; 0 1 0; 0 0 1])
    Je = sparse([1 1 0; 1 1 1; 0 1 1])
    Gp = sparse([0 6 0; 6 0 0; 0 0 0])

    @test MotifCluster.build_G(adj_mat) == G
    @test MotifCluster.build_J(adj_mat) == J
    @test MotifCluster.build_Gs(adj_mat) == Gs
    @test MotifCluster.build_Js(adj_mat) == Js
    @test MotifCluster.build_Gd(adj_mat) == Gd
    @test MotifCluster.build_Jd(adj_mat) == Jd
    @test MotifCluster.build_J0(adj_mat) == J0
    @test MotifCluster.build_Jn(adj_mat) == Jn
    @test MotifCluster.build_Id(adj_mat) == Id
    @test MotifCluster.build_Je(adj_mat) == Je
    @test MotifCluster.build_Gp(adj_mat) == Gp
end
