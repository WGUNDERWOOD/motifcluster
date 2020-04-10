#' Compute a right-multiplication with the ones matrix
#'
#' Compute \code{a * (b \%*\% one_mat)} where \code{a}, \code{b},
#' \code{ones_mat} are square matrices of the same size,
#' and \code{ones_mat} contains all entries equal to one.
#' The product \code{*} is an entry-wise (Hadamard) product,
#' while \code{\%*\%} represents matrix multiplication.
#' This method is more efficient than the naive approach
#' when \code{a} or \code{b} are sparse.
#' @param a,b Square matrices.
#' @return The square matrix \code{a * (b \%*\% one_mat)}.
#' @importFrom Matrix Diagonal drop0
#' @keywords internal

a_b_one <- function(a, b) {

  n <- nrow(a)
  ones_vec <- rep(1, n)
  ans <- Diagonal(n, as.numeric(b %*% ones_vec)) %*% a

  return(drop0(ans))
}

#' Compute a left-multiplication with the ones matrix
#'
#' Compute \code{a * (one_mat \%*\% b)} where \code{a}, \code{b},
#' \code{ones_mat} are square matrices of the same size,
#' and \code{ones_mat} contains all entries equal to one.
#' The product \code{*} is an entry-wise (Hadamard) product,
#' while \code{\%*\%} represents matrix multiplication.
#' This method is more efficient than the naive approach
#' when \code{a} or \code{b} are sparse.
#' @param a,b Square matrices.
#' @return The square matrix \code{a * (one_mat \%*\% b)}.
#' @importFrom Matrix Diagonal drop0
#' @keywords internal

a_one_b <- function(a, b) {

  n <- nrow(a)
  ones_vec <- rep(1, n)
  ans <- a %*% Diagonal(n, as.numeric(ones_vec %*% b))

  return(drop0(ans))
}

#' Set diagonal entries to zero and sparsify
#'
#' Set the diagonal entries of a matrix to zero
#' and convert it to sparse form.
#' @param some_mat A square matrix.
#' @return A sparse-form copy of the matrix with its
#' diagonal entries set to zero.
#' @importFrom Matrix drop0
#' @keywords internal

drop0_killdiag <- function(some_mat) {

  ans <- some_mat
  diag(ans) <- 0
  ans <- drop0(ans)

  return(ans)
}

#' Get largest connected component
#'
#' Get the indices of the vertices in the largest connected
#' component of a graph from its adjacency matrix.
#' @param adj_mat An adjacency matrix of a graph.
#' @return A vector of indices corresponding to the vertices in the largest
#' connected component.
#' @importFrom igraph components graph_from_adjacency_matrix
#' @export
# TODO examples

get_largest_component <- function(adj_mat) {

  # TODO remove igraph dependency

  n <- nrow(adj_mat)
  gr <- graph_from_adjacency_matrix(adj_mat + t(adj_mat))
  comps <- components(gr)
  verts_to_keep <- (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)
}

#' Get common motif names
#'
#' Get the names of some common motifs as strings.
#' @return A vector of names (strings) of common motifs.
#' @export

get_motif_names <- function() {

  motif_names <- c("Ms", "Md")

  for (i in 1:13) {
    motif_name <- paste("M", i, sep = "")
    motif_names <- c(motif_names, motif_name)
  }

  motif_names <- c(motif_names, c("Mcoll", "Mexpa"))

  return(motif_names)
}
