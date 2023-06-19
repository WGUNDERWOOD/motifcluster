@testset verbose = true "Clustering" begin

    @testset verbose = true "sample_dsbm_unweighted" begin
        Random.seed!(9349)
        sample_weight_type = "unweighted"
        block_sizes = [2, 3]
        connection_matrix = [0.4 0.5; 0.6 0.7]
        n_reps = 2000
        weight_matrix = nothing
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                     weight_matrix, sample_weight_type)
        G_vals = [0, 1]
        @test all(g in G_vals for g in G)
        G = zeros(5, 5)
        for rep in 1:n_reps
            G += MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                          weight_matrix, sample_weight_type)
        end
        G /= n_reps
        ans = sparse([0   0.4 0.5 0.5 0.5;
                      0.4   0 0.5 0.5 0.5;
                      0.6 0.6   0 0.7 0.7;
                      0.6 0.6 0.7   0 0.7;
                      0.6 0.6 0.7 0.7   0])
        @test isapprox(G, ans, atol = 0.05)
    end

    @testset verbose = true "sample_dsbm_constant_weighted" begin
        Random.seed!(2839)
        sample_weight_type = "constant"
        block_sizes = [2, 3]
        connection_matrix = [0.4 0.5; 0.6 0.7]
        n_reps = 2000
        weight_matrix = [20 30; 40 50]
        G = MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                     weight_matrix,sample_weight_type)
        G_vals = [0, 20, 30, 40, 50]
        @test all(g in G_vals for g in G)

        G = zeros(5, 5)
        for rep in 1:n_reps
            G += MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                          weight_matrix, sample_weight_type)
        end
        G /= n_reps
        ans = sparse([ 0  8 15 15 15;
                       8  0 15 15 15;
                      24 24  0 35 35;
                      24 24 35  0 35;
                      24 24 35 35  0])
        @test isapprox(G, ans, atol = 3)
  end

  @testset verbose = true "sample_dsbm_poisson_weighted" begin
      Random.seed!(2838)
      sample_weight_type = "poisson"
      block_sizes = [2, 3]
      connection_matrix = [0.4 0.5; 0.6 0.7]
      n_reps = 2000
      weight_matrix = [20 30; 40 50]
      G = MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                   weight_matrix, sample_weight_type)
      @test all(G .== floor.(G))
      @test all(G .>= 0)
      G = zeros(5, 5)
      for rep in 1:n_reps
          G += MotifCluster.sample_dsbm(block_sizes, connection_matrix,
                                        weight_matrix, sample_weight_type)
      end
      G /= n_reps
      ans = sparse([ 0  8 15 15 15;
                     8  0 15 15 15;
                    24 24  0 35 35;
                    24 24 35  0 35;
                    24 24 35 35  0])
      @test isapprox(G, ans, atol = 3)
  end

  #=

@testset verbose = true "sample_dsbm_large" begin

  rd.seed(seed = 2238)
  random.seed(2238)

  n = int(1e5)

  block_sizes = [n]
  connection_matrix = np.array([10 / n]).reshape((1, 1))

  G = MotifCluster.sample_dsbm(block_sizes, connection_matrix)

  assert G.shape == (n, n)

  end

@testset verbose = true "sample_bsbm_unweighted" begin

  rd.seed(seed = 9423)
  random.seed(9423)

  sample_weight_type = "unweighted"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 300
  bipartite_weight_matrix = None

  G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G_vals = [0, 1]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((6, 6))

  for rep in 1:n_reps
    G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G = np.array(G) / n_reps

  ans = np.array([
                         0, 0, 0, 0.3, 0.4, 0.5,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 0.05)

  return

  end


@testset verbose = true "sample_bsbm_constant_weighted" begin

  rd.seed(seed = 7482)
  random.seed(7482)

  sample_weight_type = "constant"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 300
  bipartite_weight_matrix = np.array([10, 20, 30, 40, 50, 60]).reshape((2, 3))

  G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G_vals = [0, 10, 20, 30, 40, 50, 60]
  assert (np.isin(G.toarray(), G_vals)).all()

  G = np.zeros((6, 6))

  for rep in 1:n_reps
    G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G = np.array(G) / n_reps

  ans = np.array([
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 3)

  return

  end

@testset verbose = true "sample_bsbm_poisson_weighted" begin


  rd.seed(seed = 7482)
  random.seed(7482)

  sample_weight_type = "poisson"
  source_block_sizes = [1, 2]
  dest_block_sizes = [1, 1, 1]
  bipartite_connection_matrix = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]).reshape((2, 3))
  n_reps = 300
  bipartite_weight_matrix = np.array([10, 20, 30, 40, 50, 60]).reshape((2, 3))

  G = MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  assert (G.toarray() == np.floor(G.toarray())).all()
  assert (G.toarray() >= 0).all()

  G = np.zeros((6, 6))

  for rep in 1:n_reps
    G += MotifCluster.sample_bsbm(source_block_sizes, dest_block_sizes, bipartite_connection_matrix,
                       bipartite_weight_matrix, sample_weight_type)

  G = np.array(G) / n_reps

  ans = np.array([
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
       ]).reshape((6, 6))

  assert np.allclose(G, ans, atol = 3)

  return
  end

  =#

  end
