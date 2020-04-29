context("Clustering")
library(mclust)

test_that("cluster_spectrum returns correct clusters", {

  set.seed(2352)
  n <- 10
  vects_1 <- cbind(rnorm(n, 0, 1),
                   rnorm(n, 0, 1),
                   rnorm(n, 0, 1))
  vects_2 <- cbind(rnorm(n, 0, 1),
                   rnorm(n, 4, 1),
                   rnorm(n, 4, 1))
  spectrum <- list(vects = rbind(vects_1, vects_2))

  clust_ans <- c(rep(1, n), rep(2, n))
  clust <- cluster_spectrum(spectrum, 2)

  if (clust[1] == 2) {
    clust <- 3 - clust
  }

  expect_equal(clust, clust_ans)
})

test_that("run_motif_clustering returns correct values", {

  set.seed(3957)
  n <- 50
  block_sizes <- rep(n, 3)

  connection_matrix <- matrix(c(
    0.9, 0.4, 0.4,
    0.4, 0.9, 0.4,
    0.4, 0.4, 0.9
  ), nrow = 3)

  weight_matrix <- matrix(c(
    9, 3, 3,
    3, 9, 3,
    3, 3, 9
  ), nrow = 3)

  motif_type <- "func"
  num_eigs <- 3
  num_clusts <- 3

  for (sample_weight_type in c("unweighted", "constant", "poisson")) {
    for (motif_name in get_motif_names()[1:15]) {
      for (mam_weight_type in c("unweighted", "mean", "product")) {
        for (type_lap in c("comb", "rw")) {

          # sample a new graph
          adj_mat <- sample_dsbm(block_sizes, connection_matrix,
                      weight_matrix, sample_weight_type)

          # run full method
          motif_clust_list <- run_motif_clustering(adj_mat, motif_name,
                               motif_type, mam_weight_type, "dense", num_eigs,
                               type_lap, TRUE, num_clusts)

          clusts <- motif_clust_list$clusts
          comps <- motif_clust_list$comps

          # answers
          ans_clusts <- c(rep(1, n),
                          rep(2, n),
                          rep(3, n))
          ans_clusts <- ans_clusts[comps]

          # score
          ari_score <- adjustedRandIndex(clusts, ans_clusts)

          expect_equal(ari_score, 1)
        }
      }
    }
  }
})
