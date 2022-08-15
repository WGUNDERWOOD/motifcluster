context("Sampling")

# sample_dsbm

test_that("sample_dsbm returns correct unweighted adjacency matrix", {

  set.seed(9689)

  sample_weight_type <- "unweighted"
  block_sizes <- c(2, 3)
  connection_matrix <- matrix(c(0.4, 0.5, 0.6, 0.7), nrow = 2, byrow = TRUE)
  n_reps <- 200
  weight_matrix <- NULL

  G <- sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)
  expect_true(all((G == 0) | (G == 1)))

  G <- matrix(0, nrow = 5, ncol = 5)

  for (rep in 1:n_reps) {
    G <- G + sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                         sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                           0, 0.4, 0.5, 0.5, 0.5,
                         0.4,   0, 0.5, 0.5, 0.5,
                         0.6, 0.6,   0, 0.7, 0.7,
                         0.6, 0.6, 0.7,   0, 0.7,
                         0.6, 0.6, 0.7, 0.7,   0
                       ), nrow = 5, byrow = TRUE))

  expect_equal(G, ans, tolerance = 0.05)
})

test_that("sample_dsbm returns correct constant
           weighted adjacency matrix", {

  set.seed(2239)

  sample_weight_type <- "constant"
  block_sizes <- c(2, 3)
  connection_matrix <- matrix(c(0.4, 0.5, 0.6, 0.7), nrow = 2, byrow = TRUE)
  n_reps <- 200
  weight_matrix <- matrix(c(20, 30, 40, 50), nrow = 2, byrow = TRUE)

  G <- sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)
  expect_true(all((G == 0) | (G == 20) | (G == 30) | (G == 40) | (G == 50)))

  G <- matrix(0, nrow = 5, ncol = 5)

  for (rep in 1:n_reps) {
    G <- G + sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                         sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                          0,  8, 15, 15, 15,
                          8,  0, 15, 15, 15,
                         24, 24,  0, 35, 35,
                         24, 24, 35,  0, 35,
                         24, 24, 35, 35,  0
                       ), nrow = 5, byrow = TRUE))

  expect_equal(G, ans, tolerance = 0.05)
})

test_that("sample_dsbm returns correct poisson weighted adjacency matrix", {

  set.seed(9323)

  sample_weight_type <- "poisson"
  block_sizes <- c(2, 3)
  connection_matrix <- matrix(c(0.4, 0.5, 0.6, 0.7), nrow = 2, byrow = TRUE)
  n_reps <- 200
  weight_matrix <- matrix(c(20, 30, 40, 50), nrow = 2, byrow = TRUE)

  G <- sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                   sample_weight_type)
  expect_true(all((G == floor(G)) & (G >= 0)))

  G <- matrix(0, nrow = 5, ncol = 5)

  for (rep in 1:n_reps) {
    G <- G + sample_dsbm(block_sizes, connection_matrix, weight_matrix,
                         sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                          0,  8, 15, 15, 15,
                          8,  0, 15, 15, 15,
                         24, 24,  0, 35, 35,
                         24, 24, 35,  0, 35,
                         24, 24, 35, 35,  0
                       ), nrow = 5, byrow = TRUE))

  expect_equal(G, ans, tolerance = 0.05)
})

test_that("sample_dsbm works with large sparse matrix", {

  set.seed(2238)

  n <- 1e4

  block_sizes <- c(n)
  connection_matrix <- matrix(10 / n)

  G <- sample_dsbm(block_sizes, connection_matrix)

  expect_equal(nrow(G), n)
  expect_equal(ncol(G), n)
})


# sample_bsbm

test_that("sample_bsbm returns correct unweighted adjacency matrix", {

  set.seed(9423)

  sample_weight_type <- "unweighted"
  source_block_sizes <- c(1, 2)
  dest_block_sizes <- c(1, 1, 1)
  bipartite_connection_matrix <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                        nrow = 2, byrow = TRUE)
  n_reps <- 200
  bipartite_weight_matrix <- NULL

  G <- sample_bsbm(source_block_sizes, dest_block_sizes,
                  bipartite_connection_matrix, bipartite_weight_matrix,
                  sample_weight_type)
  expect_true(all((G == 0) | (G == 1)))

  G <- matrix(0, nrow = 6, ncol = 6)

  for (rep in 1:n_reps) {
    G <- G + sample_bsbm(source_block_sizes, dest_block_sizes,
                        bipartite_connection_matrix, bipartite_weight_matrix,
                        sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                         0, 0, 0, 0.3, 0.4, 0.5,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0, 0.6, 0.7, 0.8,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0,
                         0, 0, 0,   0,   0,   0
                       ), nrow = 6, byrow = TRUE))

  expect_true(all(abs(G - ans) <= 0.1 * ans))
})

test_that("sample_bsbm returns correct constant weighted
           adjacency matrix", {

  set.seed(7482)

  sample_weight_type <- "constant"
  source_block_sizes <- c(1, 2)
  dest_block_sizes <- c(1, 1, 1)
  bipartite_connection_matrix <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                        nrow = 2, byrow = TRUE)
  n_reps <- 200
  bipartite_weight_matrix <- matrix(c(10, 20, 30, 40, 50, 60),
                                    nrow = 2, byrow = TRUE)

  G <- sample_bsbm(source_block_sizes, dest_block_sizes,
                  bipartite_connection_matrix, bipartite_weight_matrix,
                  sample_weight_type)
  expect_true(all((G == 0) | (G == 10) | (G == 20) | (G == 30)
                  | (G == 40) | (G == 50) | (G == 60)))

  G <- matrix(0, nrow = 6, ncol = 6)

  for (rep in 1:n_reps) {
    G <- G + sample_bsbm(source_block_sizes, dest_block_sizes,
                         bipartite_connection_matrix, bipartite_weight_matrix,
                         sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
                       ), nrow = 6, byrow = TRUE))

  expect_true(all(abs(G - ans) <= 0.1 * ans))
})

test_that("sample_bsbm returns correct poisson weighted adjacency matrix", {

  set.seed(7482)

  sample_weight_type <- "poisson"
  source_block_sizes <- c(1, 2)
  dest_block_sizes <- c(1, 1, 1)
  bipartite_connection_matrix <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                        nrow = 2, byrow = TRUE)
  n_reps <- 400
  bipartite_weight_matrix <- matrix(c(10, 20, 30, 40, 50, 60),
                                    nrow = 2, byrow = TRUE)

  G <- sample_bsbm(source_block_sizes, dest_block_sizes,
                  bipartite_connection_matrix, bipartite_weight_matrix,
                  sample_weight_type)
  expect_true(all((G == floor(G)) & (G >= 0)))

  G <- matrix(0, nrow = 6, ncol = 6)

  for (rep in 1:n_reps) {
    G <- G + sample_bsbm(source_block_sizes, dest_block_sizes,
                         bipartite_connection_matrix, bipartite_weight_matrix,
                         sample_weight_type) / n_reps
  }

  G <- drop0(G)
  ans <- drop0(matrix(c(
                         0, 0, 0,  3,  8, 15,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0, 24, 35, 48,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0,
                         0, 0, 0,  0,  0,  0
                       ), nrow = 6, byrow = TRUE))

  expect_true(all(abs(G - ans) <= 0.1 * ans))
})
