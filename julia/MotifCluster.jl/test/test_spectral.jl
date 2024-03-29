@testset verbose = true "Spectral" begin
    @testset verbose = true "get_first_eigs" begin
        G = sparse([7 -4 14 0; -4 19 10 0; 14 10 10 0; 0 0 0 100]')
        ans_vals = [-9, 18, 27]
        ans_vects = [-2 -1 2 0; -2 2 -1 0; -1 -2 -2 0]' ./ 3
        spect = MotifCluster.get_first_eigs(G, 3)
        vects = spect["vects"]
        vals = spect["vals"]
        for i in 1:length(vals)
            if sign(vects[1, i]) != sign(ans_vects[1, i])
                vects[:, i] = -vects[:, i]
            end
        end
        @test isapprox(vals, ans_vals)
        @test isapprox(vects, ans_vects)
    end

    @testset verbose = true "build_laplacian" begin
        G_dense = reshape(collect(0:8), (3, 3))
        G_sparse = sparse(G_dense)
        for G in [G_sparse, G_dense]
            G = G .+ G'
            degs_mat = diagm([12, 24, 36])
            comb_lap = degs_mat .- G
            rw_lap = inv(degs_mat) * (degs_mat .- G)
            @test isapprox(MotifCluster.build_laplacian(G, "comb"), comb_lap)
            @test isapprox(MotifCluster.build_laplacian(G, "rw"), rw_lap)
            G = sparse([0 1; 0 2])
            @test_throws AssertionError MotifCluster.build_laplacian(G, "rw")
        end
    end

    @testset verbose = true "run_laplace_embedding" begin
        Random.seed!(9235)
        G_dense = reshape(collect(0:8), (3, 3))
        G_sparse = sparse(G_dense)
        for G in [G_sparse, G_dense]
            G = G .+ G'
            ans_vals_comb = [0, 17.07]
            ans_vects_comb = [0.577 0.789; 0.577 -0.577; 0.577 -0.211]
            ans_vals_rw = [0, 1]
            ans_vects_rw = [0.577 0.408; 0.577 -0.816; 0.577 0.408]
            spectrum_comb = MotifCluster.run_laplace_embedding(G, 2, "comb")
            spectrum_rw = MotifCluster.run_laplace_embedding(G, 2, "rw")
            vals_comb = spectrum_comb["vals"]
            vects_comb = spectrum_comb["vects"]
            vals_rw = spectrum_rw["vals"]
            vects_rw = spectrum_rw["vects"]
            for i in 1:length(vals_comb)
                if sign(vects_comb[1, i]) != sign(ans_vects_comb[1, i])
                    vects_comb[:, i] = -vects_comb[:, i]
                end
                if sign(vects_rw[1, i]) != sign(ans_vects_rw[1, i])
                    vects_rw[:, i] = -vects_rw[:, i]
                end
            end
            @test isapprox(vals_comb, ans_vals_comb, atol=0.01)
            @test isapprox(vects_comb, ans_vects_comb, atol=0.01)
            @test isapprox(vals_rw, ans_vals_rw, atol=0.01)
            @test isapprox(vects_rw, ans_vects_rw, atol=0.01)
        end
    end

    @testset verbose = true "run_mot_embedding_restrict" begin
        Random.seed!(9235)
        adj_mat_dense = [0 2 0 0; 0 0 3 0; 4 0 0 0; 0 0 0 0]
        adj_mat_sparse = sparse(adj_mat_dense)
        for adj_mat in [adj_mat_sparse, adj_mat_dense]
            ans_adj_mat = sparse([0 2 0 0; 0 0 3 0; 4 0 0 0; 0 0 0 0])
            ans_motif_adj_mat = sparse([0 2 4 0; 2 0 3 0; 4 3 0 0; 0 0 0 0])
            ans_comps = [1, 2, 3]
            ans_adj_mat_comps = sparse([0 2 0; 0 0 3; 4 0 0])
            ans_motif_adj_mat_comps = sparse([0 2 4; 2 0 3; 4 3 0])
            ans_vals = [0, 1.354]
            ans_vects = sparse([0.577 -0.544; 0.577 0.830; 0.577 -0.126])
            embedding = MotifCluster.run_motif_embedding(adj_mat, "Ms", "func", "mean", 2, "rw",
                                                         true)
            for i in 1:length(ans_vals)
                if sign(embedding["vects"][1, i]) != sign(ans_vects[1, i])
                    embedding["vects"][:, i] = -embedding["vects"][:, i]
                end
            end
            @test isapprox(ans_adj_mat, embedding["adj_mat"])
            @test isapprox(ans_motif_adj_mat, embedding["motif_adj_mat"])
            @test isapprox(ans_comps, embedding["comps"])
            @test isapprox(ans_adj_mat_comps, embedding["adj_mat_comps"])
            @test isapprox(ans_motif_adj_mat_comps, embedding["motif_adj_mat_comps"])
            @test isapprox(ans_vals, embedding["vals"], atol=0.01)
            @test isapprox(ans_vects, embedding["vects"], atol=0.01)
        end
    end

    @testset verbose = true "run_mot_embedding_no_restrict" begin
        Random.seed!(9235)
        adj_mat_dense = [0 2 0; 0 0 3; 4 0 0]
        adj_mat_sparse = sparse(adj_mat_dense)
        for adj_mat in [adj_mat_sparse, adj_mat_dense]
            ans_adj_mat = sparse([0 2 0; 0 0 3; 4 0 0])
            ans_motif_adj_mat = sparse([0 2 4; 2 0 3; 4 3 0])
            ans_vals = [0, 1.354]
            ans_vects = sparse([0.577 -0.544; 0.577 0.830; 0.577 -0.126])
            embedding = MotifCluster.run_motif_embedding(adj_mat, "Ms", "func", "mean", 2, "rw",
                                                         false)
            for i in 1:length(ans_vals)
                if sign(embedding["vects"][1, i]) != sign(ans_vects[1, i])
                    embedding["vects"][:, i] = -embedding["vects"][:, i]
                end
            end
            @test isapprox(ans_adj_mat, embedding["adj_mat"])
            @test isapprox(ans_motif_adj_mat, embedding["motif_adj_mat"])
            @test isapprox(ans_vals, embedding["vals"], atol=0.01)
            @test isapprox(ans_vects, embedding["vects"], atol=0.01)
        end
    end
end
