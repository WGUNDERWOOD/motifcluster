@testset verbose = true "Clustering" begin

    @testset verbose = true "cluster_spectrum" begin

        Random.seed!(1)
        n = 10
        vects_1 = rand(Normal(), (n, 3))
        vects_2 = hcat(rand(Normal(), (n, 1)), rand(Normal(4), (n, 2)))

        vects =  vcat(vects_1, vects_2)
        spectrum = MotifCluster.Spectrum(vects)

        clust_ans = [[1 for _ in 1:n]; [2 for _ in 1:n]]
        clust = MotifCluster.cluster_spectrum(spectrum, 2)

        if clust[1] == 2
            clust .= 3 .- clust
        end

        println(clust)
        println(clust_ans)

        @test clust == clust_ans
    end

end


#=
def test_run_motif_clustering():

  rd.seed(3957)
  random.seed(3957)

  n = 50
  block_sizes = 3 * [n]

  connection_matrix = np.array([
    0.9, 0.4, 0.4,
    0.4, 0.9, 0.4,
    0.4, 0.4, 0.9
  ]).reshape((3, 3))

  weight_matrix = np.array([
    9, 3, 3,
    3, 9, 3,
    3, 3, 9
  ]).reshape((3, 3))

  motif_type = "func"
  num_eigs = 3
  num_clusts = 3

  for sample_weight_type in ["unweighted", "constant", "poisson"]:
    for motif_name in mcut.get_motif_names()[0:15]:
      for mam_weight_type in ["unweighted", "mean", "product"]:
        for type_lap in ["comb", "rw"]:

          # sample a new graph
          adj_mat = mcsa.sample_dsbm(block_sizes, connection_matrix,
                                     weight_matrix, sample_weight_type)

          # run full method
          motif_clust_list = mccl.run_motif_clustering(adj_mat, motif_name,
                                                       motif_type, mam_weight_type,
                                                       "dense", num_eigs,
                                                       type_lap, num_clusts,
                                                       gr_method="dense")

          clusts = motif_clust_list["clusts"]
          comps = motif_clust_list["comps"]

          # answers
          ans_clusts = n * [0] + n * [1] + n * [2]
          ans_clusts = np.take(ans_clusts, comps)

          # score
          ari_score = adjusted_rand_score(clusts, ans_clusts)

          assert ari_score == 1
=#

