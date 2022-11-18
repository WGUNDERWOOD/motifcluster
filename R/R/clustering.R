#' kmeans++ clustering
#'
#' Use the kmeans++ algorithm to cluster points
#' into \code{k} clusters, as implemented in the
#' deprecated LICORS package, using the
#' built-in function \link[stats]{kmeans}.
#' @param data An \eqn{N \times d} matrix, where there are \eqn{N} samples
#' in dimension \eqn{d}.
#' @param k The number of clusters.
#' @param iter.max The maximum number of iterations.
#' @param nstart The number of restarts.
#' @param ... Additional arguments passed to \code{\link[stats]{kmeans}}.
#' @return A list with 9 entries:
#' \itemize{
#'   \item \code{cluster}: A vector of integers from 1:k indicating the
#'     cluster to which each point is allocated.
#'   \item \code{centers}: A matrix of cluster centers.
#'   \item \code{totss}: The total sum of squares.
#'   \item \code{withinss}: Vector of within-cluster sum of squares,
#'     one component per cluster.
#'   \item \code{tot.withinss}: Total within-cluster sum of squares,
#'     i.e.sum(withinss).
#'   \item \code{betweenss}: The between-cluster sum of squares,
#'     i.e.totss-tot.withinss.
#'   \item \code{size}: The number of points in each cluster.
#'   \item \code{iter}: The number of (outer) iterations.
#'   \item \code{ifault}: An integer indicator of a possible algorithm problem.
#'   \item \code{initial.centers}: The initial centers used.
#' }
#' @importFrom stats kmeans cov
#' @examples
#' set.seed(1984)
#' n <- 100
#' X = matrix(rnorm(n), ncol = 2)
#' Y = matrix(runif(length(X)*2, -1, 1), ncol = ncol(X))
#' Z = rbind(X, Y)
#' cluster_Z = kmeanspp(Z, k = 5)
#' @seealso \code{\link[stats]{kmeans}}
#' @references
#' Arthur, D. and S. Vassilvitskii (2007).
#' ``k-means++: The advantages of careful seeding.''
#' In H. Gabow (Ed.), Proceedings of the 18th Annual ACM-SIAM
#' Symposium on Discrete Algorithms
#' [SODA07], Philadelphia, pp. 1027-1035.
#' Society for Industrial and Applied Mathematics.
#' @export

kmeanspp <- function(data, k = 2, iter.max = 100, nstart = 10, ...) {

  kk <- k

  if (length(dim(data)) == 0) {
    data <- matrix(data, ncol = 1)
  } else {
    data <- cbind(data)
  }

  num.samples <- nrow(data)
  ndim <- ncol(data)

  data.avg <- colMeans(data)
  data.cov <- cov(data)

  out <- list()
  out$tot.withinss <- Inf
  for (restart in seq_len(nstart)) {
    center_ids <- rep(0, length = kk)
    center_ids[1:2] <- sample.int(num.samples, 1)
    for (ii in 2:kk) { # the plus-plus step in kmeans
      if (ndim == 1) {
        dists <- apply(cbind(data[center_ids, ]), 1,
                       function(center) {
                         rowSums((data - rep(center, each = nrow(data)))^2)
                       })
      } else {
        dists <- apply(data[center_ids, ], 1,
                       function(center) {
                         rowSums((data - rep(center, each = nrow(data)))^2)
                       })
      }
      probs <- apply(dists, 1, min)
      probs[center_ids] <- 0
      center_ids[ii] <- sample.int(num.samples, 1, prob = probs)
    }

    tmp.out <- kmeans(data, centers = data[center_ids, ],
                      iter.max = iter.max, ...)
    tmp.out$initial.centers <- data[center_ids, ]
    if (tmp.out$tot.withinss < out$tot.withinss) {
      out <- tmp.out
    }
  }
  invisible(out)
}

#' Get cluster assignments from spectrum using k-means++
#'
#' Get a vector of cluster assignments from a spectrum,
#' using k-means++ and \code{num_clusts} clusters.
#' @param spectrum A list containing \code{$vects};
#' the matrix of eigenvectors to pass to k-means++.
#' @param num_clusts The number of clusters to find.
#' @return A length-nrow(spectrum$vects) vector of integers
#' from 1 to num_clusts, representing cluster assignments.
#' @keywords internal

cluster_spectrum <- function(spectrum, num_clusts) {

  vects <- spectrum$vects[, -1, drop = FALSE]
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
#' adj_mat <- matrix(c(1:16), nrow = 4)
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
