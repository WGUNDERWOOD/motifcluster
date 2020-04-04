context("Sampling")

# sample_dsbm

test_that("sample_dsbm returns correct unweighted adjacency matrix", {

  set.seed(9687)

  weight_type = "unweighted"
  block_sizes = c(2,3)
  connection_matrix = matrix(c(0.4, 0.5, 0.6, 0.7), nrow=length(block_sizes))
  n_reps = 1000
  weight_matrix = NULL

  G = sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)
  expect_true(all((G==0) | (G==1)))

  n = sum(block_sizes)
  G = matrix(0, nrow=n, ncol=n)

  for(rep in 1:n_reps){
    G = G + sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)/n_reps
  }
  G = unname(Matrix::drop0(G))

  ans = unname(Matrix::drop0(matrix(c(
                         0,0.4,0.6,0.6,0.6,
                         0.4,0,0.6,0.6,0.6,
                         0.5,0.5,0,0.7,0.7,
                         0.5,0.5,0.7,0,0.7,
                         0.5,0.5,0.7,0.7,0
                       ), nrow=5, byrow=TRUE)))

  expect_equal(G, ans, tolerance=0.05)
})

test_that("sample_dsbm returns correct deterministic weighted adjacency matrix", {

  set.seed(2239)

  weight_type = "deterministic"
  block_sizes = c(2,3)
  connection_matrix = matrix(c(0.4, 0.5, 0.6, 0.7), nrow=length(block_sizes))
  n_reps = 1000
  weight_matrix = matrix(c(20, 30, 40, 50), nrow=length(block_sizes))

  G = sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)
  expect_true(all((G==0) | (G==20) | (G==30) | (G==40) | (G==50)))

  n = sum(block_sizes)
  G = matrix(0, nrow=n, ncol=n)

  for(rep in 1:n_reps){
    G = G + sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)/n_reps
  }
  G = unname(Matrix::drop0(G))

  ans = unname(Matrix::drop0(matrix(c(
                         0,8,24,24,24,
                         8,0,24,24,24,
                         15,15,0,35,35,
                         15,15,35,0,35,
                         15,15,35,35,0
                       ), nrow=5, byrow=TRUE)))

  expect_equal(G, ans, tolerance=0.05)
})

test_that("sample_dsbm returns correct poisson weighted adjacency matrix", {

  set.seed(9323)

  weight_type = "poisson"
  block_sizes = c(2,3)
  connection_matrix = matrix(c(0.4, 0.5, 0.6, 0.7), nrow=length(block_sizes))
  n_reps = 1000
  weight_matrix = matrix(c(20, 30, 40, 50), nrow=length(block_sizes))

  G = sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)
  expect_true(all((G==floor(G)) & (G>=0)))

  n = sum(block_sizes)
  G = matrix(0, nrow=n, ncol=n)

  for(rep in 1:n_reps){
    G = G + sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)/n_reps
  }
  G = unname(Matrix::drop0(G))

  ans = unname(Matrix::drop0(matrix(c(
                         0,8,24,24,24,
                         8,0,24,24,24,
                         15,15,0,35,35,
                         15,15,35,0,35,
                         15,15,35,35,0
                       ), nrow=5, byrow=TRUE)))

  expect_equal(G, ans, tolerance=0.05)
})
