@testset verbose = true "Spectral" begin

    @testset verbose = true "get_first_eigs" begin
        G = sparse([7 -4 14 0; -4 19 10 0; 14 10 10 0; 0 0 0 100]')
        ans_vals = [-9, 18, 27]
        ans_vects = [-2 -1 2 0; -2 2 -1 0; -1 -2 -2 0]' ./ 3
        spect = MotifCluster.get_first_eigs(G, 3)
        vects = spect.vects
        vals = spect.vals
        for i in 1:length(vals)
            if sign(vects[1, i]) != sign(ans_vects[1, i])
                vects[:, i] = -vects[:, i]
            end
        end
        @test isapprox(vals, ans_vals)
        @test isapprox(vects, ans_vects)
    end

    @testset verbose = true "build_laplacian" begin
        G = sparse(reshape(collect(0:8), (3, 3)))
        G = G .+ G'
        degs_mat = diagm([12, 24, 36])
        comb_lap = degs_mat .- G
        rw_lap = inv(degs_mat) * (degs_mat .- G)
        @test isapprox(MotifCluster.build_laplacian(G, "comb"), comb_lap)
        @test isapprox(MotifCluster.build_laplacian(G, "rw"), rw_lap)
        G = sparse([0 1; 0 2])
        @test_throws AssertionError MotifCluster.build_laplacian(G, "rw")
    end

    @testset verbose = true "run_laplace_embedding" begin
        Random.seed!(9235)
        G = sparse(reshape(collect(0:8), (3, 3)))
        G = G .+ G'
        ans_vals_comb = [0, 17.07]
        ans_vects_comb = [0.577 0.789; 0.577 -0.577; 0.577 -0.211]
        ans_vals_rw = [0, 1]
        ans_vects_rw = [0.577 0.408; 0.577 -0.816; 0.577 0.408]
        spectrum_comb = MotifCluster.run_laplace_embedding(G, 2, "comb")
        spectrum_rw = MotifCluster.run_laplace_embedding(G, 2, "rw")
        vals_comb = spectrum_comb.vals
        vects_comb = spectrum_comb.vects
        vals_rw = spectrum_rw.vals
        vects_rw = spectrum_rw.vects
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

    #=

    # run_motif_embedding

    def test_run_mot_embedding_dense_restrict():

    np.random.seed(9235)

    adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4))

    # answers
    ans_adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4))

    ans_motif_adj_mat = np.array([
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4))

    ans_comps = [0, 1, 2]

    ans_adj_mat_comps = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3))

    ans_motif_adj_mat_comps = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
    ]).reshape((3, 3))

    ans_vals = [0, 1.354]

    ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
    ]).reshape((3, 2))

    # run motif embedding
    emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
    "rw", restrict=True)

    # flip eigenvector signs if necessary
    for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
    emb_list["vects"][:, i] = -emb_list["vects"][:, i]

    assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
    assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
    assert np.allclose(ans_comps, emb_list["comps"])
    assert np.allclose(ans_adj_mat_comps, emb_list["adj_mat_comps"].toarray())
    assert np.allclose(ans_motif_adj_mat_comps, emb_list["motif_adj_mat_comps"].toarray())
    assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
    assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


    def test_run_mot_embedding_sparse_restrict():

    np.random.seed(9235)

    adj_mat = sparse.csr_matrix(np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4)))

    # answers
    ans_adj_mat = np.array([
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4))

    ans_motif_adj_mat = np.array([
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
    ]).reshape((4, 4))

    ans_comps = [0, 1, 2]

    ans_adj_mat_comps = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3))

    ans_motif_adj_mat_comps = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
    ]).reshape((3, 3))

    ans_vals = [0, 1.354]

    ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
    ]).reshape((3, 2))

    # run motif embedding
    emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
    "rw", restrict=True)

    # flip eigenvector signs if necessary
    for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
    emb_list["vects"][:, i] = -emb_list["vects"][:, i]

    assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
    assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
    assert np.allclose(ans_comps, emb_list["comps"])
    assert np.allclose(ans_adj_mat_comps, emb_list["adj_mat_comps"].toarray())
    assert np.allclose(ans_motif_adj_mat_comps, emb_list["motif_adj_mat_comps"].toarray())
    assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
    assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


    def test_run_mot_embedding_dense_no_restrict():

    np.random.seed(9235)

    adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3))

    # answers
    ans_adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3))

    ans_motif_adj_mat = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
    ]).reshape((3, 3))

    ans_vals = [0, 1.354]

    ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
    ]).reshape((3, 2))

    # run motif embedding
    emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
    "rw", restrict=False)

    # flip eigenvector signs if necessary
    for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
    emb_list["vects"][:, i] = -emb_list["vects"][:, i]

    assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
    assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
    assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
    assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)


    def test_run_mot_embedding_sparse_no_restrict():

    np.random.seed(9235)

    adj_mat = sparse.csr_matrix(np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3)))

    # answers
    ans_adj_mat = np.array([
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
    ]).reshape((3, 3))

    ans_motif_adj_mat = np.array([
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
    ]).reshape((3, 3))

    ans_vals = [0, 1.354]

    ans_vects = np.array([
    0.577, -0.544,
    0.577, 0.830,
    0.577, -0.126
    ]).reshape((3, 2))

    # run motif embedding
    emb_list = mcsp.run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
    "rw", restrict=False)

    # flip eigenvector signs if necessary
    for i in range(len(ans_vals)):
    if np.sign(emb_list["vects"][0, i]) != np.sign(ans_vects[0, i]):
    emb_list["vects"][:, i] = -emb_list["vects"][:, i]

    assert np.allclose(ans_adj_mat, emb_list["adj_mat"].toarray())
    assert np.allclose(ans_motif_adj_mat, emb_list["motif_adj_mat"].toarray())
    assert np.allclose(ans_vals, emb_list["vals"], atol=0.01)
    assert np.allclose(ans_vects, emb_list["vects"], atol=0.01)

    =#

end
