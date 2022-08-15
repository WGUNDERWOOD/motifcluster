context("Spectral methods")

# get_first_eigs

test_that("get_first_eigs returns correct values on dense matrix", {

  G <- matrix(c(7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100), nrow = 4)
  ans_vals <- c(-9, 18, 27)
  ans_vects <- matrix(
    c(-2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0) / 3, nrow = 4)

  spect <- get_first_eigs(G, 3)
  vals <- spect$vals
  vects <- spect$vects

  for (i in seq_len(length(vals))) {
    if (sign(vects[1, i]) != sign(ans_vects[1, i])) {
      vects[, i] <- -vects[, i]
    }
  }

  expect_equal(vals, ans_vals)
  expect_equal(vects, ans_vects)
})

test_that("get_first_eigs returns correct values on sparse matrix", {

  G <- drop0(matrix(c(7, -4, 14, 0, -4, 19, 10, 0,
                14, 10, 10, 0, 0, 0, 0, 100), nrow = 4))
  ans_vals <- c(-9, 18, 27)
  ans_vects <- matrix(c(
    -2, -1, 2, 0, -2, 2, -1, 0, -1, -2, -2, 0) / 3, nrow = 4)

  spect <- get_first_eigs(G, 3)
  vals <- spect$vals
  vects <- spect$vects

  for (i in seq_len(length(vals))) {
    if (sign(vects[1, i]) != sign(ans_vects[1, i])) {
      vects[, i] <- -vects[, i]
    }
  }

  expect_equal(vals, ans_vals)
  expect_equal(vects, ans_vects)
})

# build_laplacian

test_that("build_laplacian returns correct matrices on dense matrix", {

  G <- matrix(c(0:8), nrow = 3)
  G <- G + t(G)

  degs_mat <- diag(c(12, 24, 36))
  comb_lap <- drop0(degs_mat - G)
  rw_lap <- drop0(solve(degs_mat) %*% (degs_mat - G))

  expect_equal(build_laplacian(G, type_lap = "comb"), comb_lap)
  expect_equal(build_laplacian(G, type_lap = "rw"), rw_lap)
})

test_that("build_laplacian returns correct matrices on sparse matrix", {

  G <- matrix(c(0:8), nrow = 3)
  G <- drop0(G + t(G))

  degs_mat <- diag(c(12, 24, 36))
  comb_lap <- degs_mat - G
  rw_lap <- drop0(solve(degs_mat) %*% (degs_mat - G))

  expect_equal(build_laplacian(G, type_lap = "comb"), comb_lap)
  expect_equal(build_laplacian(G, type_lap = "rw"), rw_lap)
})

test_that("build_laplacian gives correct error if row sums are zero", {

  G <- drop0(matrix(c(0, 1, 0, 2)))
  expect_error(build_laplacian(G, type_lap = "rw"),
               "row sums of adj_mat must be non-zero")
})

# run_laplace_embedding

test_that("run_laplace_embedding returns correct spectrum on dense matrix", {

  set.seed(9235)

  G <- matrix(c(0:8), nrow = 3)
  G <- G + t(G)

  ans_vals_comb <- c(0, 17.07)
  ans_vects_comb <- matrix(c(
    0.577,  0.789,
    0.577, -0.577,
    0.577, -0.211
  ), nrow = 3, byrow = TRUE)

  ans_vals_rw <- c(0, 1)
  ans_vects_rw <- matrix(c(
    0.577,  0.408,
    0.577, -0.816,
    0.577,  0.408
  ), nrow = 3, byrow = TRUE)

  spectrum_comb <- run_laplace_embedding(G, 2, "comb")
  spectrum_rw <- run_laplace_embedding(G, 2, "rw")

  vals_comb <- spectrum_comb$vals
  vects_comb <- spectrum_comb$vects
  vals_rw <- spectrum_rw$vals
  vects_rw <- spectrum_rw$vects

  for (i in seq_len(length(vals_comb))) {
    if (sign(vects_comb[1, i]) != sign(ans_vects_comb[1, i])) {
      vects_comb[, i] <- -vects_comb[, i]
    }
    if (sign(vects_rw[1, i]) != sign(ans_vects_rw[1, i])) {
      vects_rw[, i] <- -vects_rw[, i]
    }
  }

  expect_equal(vals_comb, ans_vals_comb, tolerance = 0.01)
  expect_equal(vects_comb, ans_vects_comb, tolerance = 0.01)
  expect_equal(vals_rw, ans_vals_rw, tolerance = 0.01)
  expect_equal(vects_rw, ans_vects_rw, tolerance = 0.01)

})

test_that("run_laplace_embedding returns correct spectrum on sparse matrix", {

  set.seed(9235)

  G <- matrix(c(0:8), nrow = 3)
  G <- drop0(G + t(G))

  ans_vals_comb <- c(0, 17.07)
  ans_vects_comb <- matrix(c(
    0.577,  0.789,
    0.577, -0.577,
    0.577, -0.211
  ), nrow = 3, byrow = TRUE)

  ans_vals_rw <- c(0, 1)
  ans_vects_rw <- matrix(c(
    0.577,  0.408,
    0.577, -0.816,
    0.577,  0.408
  ), nrow = 3, byrow = TRUE)

  spectrum_comb <- run_laplace_embedding(G, 2, "comb")
  spectrum_rw <- run_laplace_embedding(G, 2, "rw")

  vals_comb <- spectrum_comb$vals
  vects_comb <- spectrum_comb$vects
  vals_rw <- spectrum_rw$vals
  vects_rw <- spectrum_rw$vects

  for (i in seq_len(length(vals_comb))) {
    if (sign(vects_comb[1, i]) != sign(ans_vects_comb[1, i])) {
      vects_comb[, i] <- -vects_comb[, i]
    }
    if (sign(vects_rw[1, i]) != sign(ans_vects_rw[1, i])) {
      vects_rw[, i] <- -vects_rw[, i]
    }
  }

  expect_equal(vals_comb, ans_vals_comb, tolerance = 0.01)
  expect_equal(vects_comb, ans_vects_comb, tolerance = 0.01)
  expect_equal(vals_rw, ans_vals_rw, tolerance = 0.01)
  expect_equal(vects_rw, ans_vects_rw, tolerance = 0.01)

})

# run_motif_embedding

test_that("run_motif_embedding correct on dense matrix with restrict", {

  set.seed(9235)

  adj_mat <- matrix(c(
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE)

  # answers
  ans_adj_mat <- drop0(matrix(c(
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE))

  ans_motif_adj_mat <- drop0(matrix(c(
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE))

  ans_comps <- 1:3

  ans_adj_mat_comps <- drop0(matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE))

  ans_motif_adj_mat_comps <- drop0(matrix(c(
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ), nrow = 3, byrow = TRUE))

  ans_vals <- c(0, 1.354)

  ans_vects <- matrix(c(
    0.577, 0.544,
    0.577, -0.830,
    0.577, 0.126
  ), nrow = 3, byrow = TRUE)

  # run motif embedding
  emb_list <- run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                 "rw", restrict = TRUE)

  # flip eigenvector signs if necessary
  for (i in seq_len(length(ans_vals))) {
    if (sign(emb_list$vects[1, i]) != sign(ans_vects[1, i])) {
      emb_list$vects[, i] <- -emb_list$vects[, i]
    }
  }

  expect_equal(ans_adj_mat, emb_list$adj_mat)
  expect_true(all(ans_motif_adj_mat - emb_list$motif_adj_mat) == 0)
  expect_equal(ans_comps, emb_list$comps)
  expect_equal(ans_adj_mat_comps, emb_list$adj_mat_comps)
  expect_true(all(ans_motif_adj_mat_comps - emb_list$motif_adj_mat_comps) == 0)
  expect_equal(ans_vals, emb_list$vals, tolerance = 0.01)
  expect_equal(ans_vects, emb_list$vects, tolerance = 0.01)
})

test_that("run_motif_embedding correct on sparse matrix with restrict", {

  set.seed(9235)

  adj_mat <- matrix(c(
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE)

  # answers
  ans_adj_mat <- drop0(matrix(c(
    0, 2, 0, 0,
    0, 0, 3, 0,
    4, 0, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE))

  ans_motif_adj_mat <- drop0(matrix(c(
    0, 2, 4, 0,
    2, 0, 3, 0,
    4, 3, 0, 0,
    0, 0, 0, 0
  ), nrow = 4, byrow = TRUE))

  ans_comps <- 1:3

  ans_adj_mat_comps <- drop0(matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE))

  ans_motif_adj_mat_comps <- drop0(matrix(c(
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ), nrow = 3, byrow = TRUE))

  ans_vals <- c(0, 1.354)

  ans_vects <- matrix(c(
    0.577, 0.544,
    0.577, -0.830,
    0.577, 0.126
  ), nrow = 3, byrow = TRUE)

  # run motif embedding
  emb_list <- run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                  "rw", restrict = TRUE)

  # flip eigenvector signs if necessary
  for (i in seq_len(length(ans_vals))) {
    if (sign(emb_list$vects[1, i]) != sign(ans_vects[1, i])) {
      emb_list$vects[, i] <- -emb_list$vects[, i]
    }
  }

  expect_equal(ans_adj_mat, emb_list$adj_mat)
  expect_true(all(ans_motif_adj_mat - emb_list$motif_adj_mat) == 0)
  expect_equal(ans_comps, emb_list$comps)
  expect_equal(ans_adj_mat_comps, emb_list$adj_mat_comps)
  expect_true(all(ans_motif_adj_mat_comps - emb_list$motif_adj_mat_comps) == 0)
  expect_equal(ans_vals, emb_list$vals, tolerance = 0.01)
  expect_equal(ans_vects, emb_list$vects, tolerance = 0.01)
})

test_that("run_motif_embedding correct on dense matrix without restrict", {

  set.seed(9235)

  adj_mat <- matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE)

  # answers
  ans_adj_mat <- drop0(matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE))

  ans_motif_adj_mat <- drop0(matrix(c(
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ), nrow = 3, byrow = TRUE))

  ans_vals <- c(0, 1.354)

  ans_vects <- matrix(c(
    0.577, 0.544,
    0.577, -0.830,
    0.577, 0.126
  ), nrow = 3, byrow = TRUE)

  # run motif embedding
  emb_list <- run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                  "rw", restrict = FALSE)

  # flip eigenvector signs if necessary
  for (i in seq_len(length(ans_vals))) {
    if (sign(emb_list$vects[1, i]) != sign(ans_vects[1, i])) {
      emb_list$vects[, i] <- -emb_list$vects[, i]
    }
  }

  expect_equal(ans_adj_mat, emb_list$adj_mat)
  expect_true(all(ans_motif_adj_mat - emb_list$motif_adj_mat) == 0)
  expect_equal(ans_vals, emb_list$vals, tolerance = 0.01)
  expect_equal(ans_vects, emb_list$vects, tolerance = 0.01)
})

test_that("run_motif_embedding correct on sparse matrix without restrict", {

  set.seed(9235)

  adj_mat <- matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE)

  # answers
  ans_adj_mat <- drop0(matrix(c(
    0, 2, 0,
    0, 0, 3,
    4, 0, 0
  ), nrow = 3, byrow = TRUE))

  ans_motif_adj_mat <- drop0(matrix(c(
    0, 2, 4,
    2, 0, 3,
    4, 3, 0
  ), nrow = 3, byrow = TRUE))

  ans_vals <- c(0, 1.354)

  ans_vects <- matrix(c(
    0.577, 0.544,
    0.577, -0.830,
    0.577, 0.126
  ), nrow = 3, byrow = TRUE)

  # run motif embedding
  emb_list <- run_motif_embedding(adj_mat, "Ms", "func", "mean", "dense", 2,
                                  "rw", restrict = FALSE)

  # flip eigenvector signs if necessary
  for (i in seq_len(length(ans_vals))) {
    if (sign(emb_list$vects[1, i]) != sign(ans_vects[1, i])) {
      emb_list$vects[, i] <- -emb_list$vects[, i]
    }
  }

  expect_equal(ans_adj_mat, emb_list$adj_mat)
  expect_true(all(ans_motif_adj_mat - emb_list$motif_adj_mat) == 0)
  expect_equal(ans_vals, emb_list$vals, tolerance = 0.01)
  expect_equal(ans_vects, emb_list$vects, tolerance = 0.01)
})
