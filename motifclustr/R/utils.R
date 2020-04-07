#' Compute a right-multiplication with the ones matrix
#'
#' Compute A * (B %*% one_mat) where A, B, ones_mat are square matrices of the same size,
#' and one_mat is a ones matrix. The product * is an entry-wise (Hadamard) product,
#' while %*% represents matrix multiplication.
#' This method is more efficient than the naive approach when A or B are sparse.
#' @param A A square n by n matrix.
#' @param B A square n by n matrix.
#' @return A * (B %*% one_mat)
#' @importFrom Matrix Diagonal drop0

a_b_one <- function(A, B){

  n = nrow(A)
  ones_vec = rep(1, n)
  ans = Diagonal(n, as.numeric(B %*% ones_vec)) %*% A

  return(drop0(ans))
}

#' Compute a left-multiplication with the ones matrix
#'
#' Compute A * (one_mat %*% B) where A, B, ones_mat are square matrices of the same size,
#' and one_mat is a ones matrix. The product * is an entry-wise (Hadamard) product,
#' while %*% represents matrix multiplication.
#' This method is more efficient than the naive approach when A or B are sparse.
#' @param A A square n by n matrix.
#' @param B A square n by n matrix.
#' @return A * (one_mat %*% B)
#' @importFrom Matrix Diagonal drop0

a_one_b <- function(A, B){

  n = nrow(A)
  ones_vec = rep(1, n)
  ans = A %*% Diagonal(n, as.numeric(ones_vec %*% B))

  return(drop0(ans))
}

#' Set diagonal entries to zero and sparsify
#'
#' Set the diagonal entries of a matrix to zero
#' and convert it to sparse form.
#' @param some_mat Some matrix.
#' @return A sparse-form copy of the matrix with its diagonal entries set to zero.
#' @importFrom Matrix drop0

drop0_killdiag <- function(some_mat){

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

get_largest_component <- function(adj_mat){

  n <- nrow(adj_mat)
  gr <- graph_from_adjacency_matrix(adj_mat + t(adj_mat))
  comps <- components(gr)
  verts_to_keep <- (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)
}

#' Get common motif names
#'
#' Get the names of some common motifs as strings.
#' @return A vector of names of common motifs.

get_motif_names <- function(){

  motif_names = c("Ms", "Md")

  for(i in 1:13){
    motif_name = paste("M", i, sep="")
    motif_names = c(motif_names, motif_name)
  }

  motif_names = c(motif_names, c("Mcoll", "Mexpa"))

  return(motif_names)
}
