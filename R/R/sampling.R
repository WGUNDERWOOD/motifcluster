#' Sample a directed stochastic block model (DSBM)
#'
#' Sample the (weighted) adjacency matrix of a (weighted) directed stochastic
#' block model (DSBM) with specified parameters.
#' @param block_sizes A vector containing the size of each block of vertices.
#' @param connection_matrix A matrix containing the block-to-block connection
#' probabilities.
#' @param sample_weight_type The type of weighting scheme.
#' One of \code{"unweighted"}, \code{"constant"} or \code{"poisson"}.
#' @param weight_matrix A matrix containing the block-to-block weight
#' parameters.
#' Unused for \code{sample_weight_type = "constant"}.
#' Defaults to \code{NULL}.
#' @return A randomly sampled (weighted) adjacency matrix of a DSBM.
#' @export
#' @importFrom stats rpois rbinom
#' @examples
#' block_sizes <- c(10, 10)
#' connection_matrix <- matrix(c(0.8, 0.1, 0.1, 0.8), nrow = 2, byrow = TRUE)
#' weight_matrix <- matrix(c(10, 3, 3, 10), nrow = 2, byrow = TRUE)
#' sample_dsbm(block_sizes, connection_matrix, weight_matrix, "poisson")

sample_dsbm <- function(block_sizes, connection_matrix,
  weight_matrix = NULL,
  sample_weight_type = c("unweighted", "constant", "poisson")) {

  # check args
  sample_weight_type <- match.arg(sample_weight_type)
  stopifnot(all.equal(block_sizes, as.integer(block_sizes)))
  stopifnot(all(block_sizes > 0))
  stopifnot(length(block_sizes) == nrow(connection_matrix))
  stopifnot(length(block_sizes) == ncol(connection_matrix))
  stopifnot(all(connection_matrix >= 0))
  stopifnot(all(connection_matrix <= 1))

  if (sample_weight_type != "unweighted") {
    stopifnot(!is.null(weight_matrix))
  }

  if (!is.null(weight_matrix)) {
    stopifnot(length(block_sizes) == nrow(weight_matrix))
    stopifnot(length(block_sizes) == ncol(weight_matrix))
    stopifnot(all(weight_matrix >= 0))
  }

  # initialize variables
  k <- length(block_sizes)
  block_list <- list()

  for (i in 1:k) {

    row_list <- list()

    for (j in 1:k) {

      # block parameters
      ni <- block_sizes[i]
      nj <- block_sizes[j]
      p <- connection_matrix[i, j]

      # generate block
      if (sample_weight_type == "unweighted") {
        w <- 1
      } else {
        w <- weight_matrix[i, j]
      }

      block <- random_sparse_matrix(ni, nj, p, sample_weight_type, w)
      row_list <- append(row_list, block)
    }

    row_block <- do.call(cbind, row_list)
    block_list <- append(block_list, row_block)
  }

  adj_mat <- do.call(rbind, block_list)
  adj_mat <- drop0_killdiag(adj_mat)

  return(adj_mat)
}

#' Sample a bipartite stochastic block model (BSBM)
#'
#' Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic
#' block model (BSBM) with specified parameters.
#' @param source_block_sizes A vector containing the size of each block
#' of source vertices.
#' @param dest_block_sizes A vector containing the size of each block
#' of destination vertices.
#' @param bipartite_connection_matrix A matrix containing the
#' source block to destination block
#' connection probabilities.
#' @param sample_weight_type The type of weighting scheme.
#' One of \code{"unweighted"}, \code{"constant"} or \code{"poisson"}.
#' @param bipartite_weight_matrix A matrix containing the
#' source block to destination block weight parameters.
#' Unused for \code{sample_weight_type = "constant"}.
#' Defaults to \code{NULL}.
#' @return A randomly sampled (weighted) adjacency matrix of a BSBM.
#' @export
#' @examples
#' source_block_sizes <- c(10, 10)
#' dest_block_sizes <- c(10, 10, 10)
#' bipartite_connection_matrix <- matrix(c(0.8, 0.5, 0.1, 0.1, 0.5, 0.8),
#'       nrow = 2, byrow = TRUE)
#' bipartite_weight_matrix = matrix(c(20, 10, 2, 2, 10, 20),
#'       nrow = 2, byrow = TRUE)
#' sample_bsbm(source_block_sizes, dest_block_sizes,
#'       bipartite_connection_matrix, bipartite_weight_matrix, "poisson")

sample_bsbm <- function(source_block_sizes, dest_block_sizes,
  bipartite_connection_matrix,
  bipartite_weight_matrix = NULL,
  sample_weight_type = c("unweighted", "constant", "poisson")) {

  # check args
  sample_weight_type <- match.arg(sample_weight_type)
  stopifnot(length(source_block_sizes) == nrow(bipartite_connection_matrix))
  stopifnot(length(dest_block_sizes) == ncol(bipartite_connection_matrix))

  if (sample_weight_type != "unweighted") {
    stopifnot(!is.null(bipartite_weight_matrix))
  }

  if (!is.null(bipartite_weight_matrix)) {
    stopifnot(length(source_block_sizes) == nrow(bipartite_weight_matrix))
    stopifnot(length(dest_block_sizes) == ncol(bipartite_weight_matrix))
    stopifnot(all(bipartite_weight_matrix >= 0))
  }

  # initialize parameters
  ks <- length(source_block_sizes)
  kd <- length(dest_block_sizes)
  zeros_ss <- matrix(0, nrow = ks, ncol = ks)
  zeros_d <- matrix(0, nrow = kd, ncol = (ks + kd))

  # build block sizes vector
  block_sizes <- c(source_block_sizes, dest_block_sizes)

  # build connection matrix
  connection_matrix <- cbind(zeros_ss, bipartite_connection_matrix)
  connection_matrix <- rbind(connection_matrix, zeros_d)

  # build weight matrix
  if (!is.null(bipartite_weight_matrix)) {
    weight_matrix <- cbind(zeros_ss, bipartite_weight_matrix)
    weight_matrix <- rbind(weight_matrix, zeros_d)
  } else {
    weight_matrix <- NULL
  }

  # sample BSBM
  adj_mat <- sample_dsbm(block_sizes, connection_matrix,
                         weight_matrix, sample_weight_type)

  return(adj_mat)
}

#' Generate a small graph for demonstrations
#'
#' Generate the sparse and dense adjacency matrices of a small weighted
#' directed graph, for demonstrating methods and running tests.
#' @return A list with two entries:
#' \code{adj_mat_dense} is the adjacency matrix in dense form, and
#' \code{adj_mat_sparse} is the adjacency matrix in sparse form.
#' @keywords internal

demonstration_graph <- function() {

  adj_mat_dense <- matrix(c(
    0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,
    2, 0,  3,  0,  6,  8,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0, 10,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,
    4, 5,  0,  0,  0, 14,  0,  0, 18, 19,  0, 0,
    0, 7,  9,  0, 13,  0,  0,  0,  0,  0, 21, 0,
    0, 0, 11, 12,  0, 15,  0, 17,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 16,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0,  0,  0,  0, 24,  0, 0,
    0, 0,  0,  0,  0, 20,  0,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 22,  0,  0,  0,  0, 0,
    0, 0,  0,  0,  0,  0, 23,  0,  0,  0,  0, 0
  ), nrow = 12, byrow = TRUE)

  adj_mat_sparse <- drop0(adj_mat_dense)

  ans <- list(adj_mat_dense = adj_mat_dense, adj_mat_sparse = adj_mat_sparse)

  return(ans)
}
