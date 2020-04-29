#' Get cluster assignments from spectrum using k-means++
#'
#' Get a vector of cluster assignments from a spectrum,
#' using k-means++ and \code{num_clusts} clusters.
#' @param spectrum A list containing \code{$vects};
#' the matrix of eigenvectors to pass to k-means++.
#' @param num_clusts The number of clusters to find.
#' @return A length-nrow(spectrum$vects) vector of integers
#' from 1 to num_clusts, representing cluster assignments.
#' @importFrom LICORS kmeanspp
#' @keywords internal

cluster_spectrum <- function(spectrum, num_clusts) {

  vects <- spectrum$vects[, -1]
  kmeans_plus_plus <- kmeanspp(vects, k = num_clusts)
  cluster_assigns <- kmeans_plus_plus$cluster

  return(cluster_assigns)
}

#' Run motif-based clustering
#'
#' Run motif-based clustering on the adjacency matrix of a
#' (weighted directed) network,
#' using a specified motif, motif type, weighting scheme,
#' embedding dimension, number of clusters and Laplacian type.
#' @param adj_mat Adjacency matrix to be embedded.
#' @param motif_name Motif used for the motif adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to use.
#' One of \code{"func"} or \code{"struc"}.
#' @param mam_weight_type Weighting scheme for the motif adjacency matrix.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method The method to use for building the motif adjacency matrix.
#' One of \code{"sparse"} or \code{"dense"}.
#' @param num_eigs Number of eigenvalues and eigenvectors for the embedding.
#' @param type_lap Type of Laplacian for the embedding.
#' One of \code{"comb"} or \code{"rw"}.
#' @param restrict Whether or not to restrict the motif adjacency matrix
#' to its largest connected component before embedding.
#' @param num_clusts The number of clusters to find.
#' @return A list with 8 entries:
#' \itemize{
#'   \item \code{adj_mat}: the original adjacency matrix.
#'   \item \code{motif_adj_mat}: the motif adjacency matrix.
#'   \item \code{comps}: the indices of the largest connected component
#'     of the motif adjacency matrix
#'     (if restrict = TRUE).
#'   \item \code{adj_mat_comps}: the original adjacency matrix restricted
#'     to the largest connected component of the motif adjacency matrix
#'     (if restrict = TRUE).
#'   \item \code{motif_adj_mat_comps}: the motif adjacency matrix restricted
#'     to its largest connected component
#'     (if restrict = TRUE).
#'   \item \code{vals}: a length-\code{num_eigs} vector containing the
#'     eigenvalues associated with the Laplace embedding
#'     of the (restricted) motif adjacency matrix.
#'   \item \code{vects}: a matrix
#'     containing the eigenvectors associated with the Laplace embedding
#'     of the (restricted) motif adjacency matrix.
#'   \item \code{clusts}: a vector containing integers representing the
#'     cluster assignment of each vertex in the (restricted) graph.
#' }
#' @examples
#' adj_mat <- matrix(c(1:9), nrow = 3)
#' run_motif_clustering(adj_mat, "M1", "func")
#' @export

run_motif_clustering <- function(adj_mat, motif_name,
  motif_type = c("struc", "func"),
  mam_weight_type = c("unweighted", "mean", "product"),
  mam_method = c("sparse", "dense"),
  num_eigs = 2,
  type_lap = c("comb", "rw"),
  restrict = TRUE,
  num_clusts = 2) {

  motif_type <- match.arg(motif_type)
  mam_weight_type <- match.arg(mam_weight_type)
  mam_method <- match.arg(mam_method)
  type_lap <- match.arg(type_lap)
  stopifnot(typeof(restrict) == typeof(TRUE))

  spectrum <- run_motif_embedding(
    adj_mat, motif_name, motif_type, mam_weight_type,
    mam_method, num_eigs, type_lap, restrict)

  cluster_assigns <- cluster_spectrum(spectrum, num_clusts)

  ans <- list()
  ans$adj_mat <- adj_mat
  ans$motif_adj_mat <- spectrum$motif_adj_mat
  ans$comps <- spectrum$comps
  ans$adj_mat_comps <- spectrum$adj_mat_comps
  ans$motif_adj_mat_comps <- spectrum$motif_adj_mat_comps
  ans$vals <- spectrum$vals
  ans$vects <- spectrum$vects
  ans$clusts <- cluster_assigns

  return(ans)
}
