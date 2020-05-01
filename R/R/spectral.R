#' Compute first few eigenvalues and eigenvectors of a matrix
#'
#' Compute the first few eigenvalues (by magnitude) and
#' associated eigenvectors of a matrix.
#' @param some_mat Matrix for which eigenvalues and eigenvectors
#' are to be calculated.
#' @param num_eigs Number of eigenvalues and eigenvectors to calculate.
#' @return A list with two entries:
#' \code{vals} contains a length-\code{num_eigs} vector of the first few
#' eigenvalues,
#' and vects contains an \code{nrow(some_mat)} by \code{num_eigs} matrix
#' of the associated eigenvectors.
#' @importFrom RSpectra eigs
#' @keywords internal

get_first_eigs <- function(some_mat, num_eigs) {

  # check args
  stopifnot(all.equal(num_eigs, as.integer(num_eigs)))
  stopifnot(num_eigs > 0)

  # get spectrum for large num_eigs
  if (num_eigs >= nrow(some_mat) - 1) {
    spectrum_eigs <- eigen(some_mat)
  }
  else {
    spectrum_eigs <- eigs(some_mat, num_eigs, which = "SM")
  }

  # use real parts
  vals <- Re(spectrum_eigs$values)
  vects <- Re(spectrum_eigs$vectors)

  # order eigenvalues and eigenvectors
  ordering <- order(vals)
  vals <- vals[ordering]
  vects <- vects[, ordering]

  # only return the specified number (eigen may return more)
  vals <- vals[1:num_eigs]
  vects <- vects[, 1:num_eigs]

  # return a list
  spectrum <- list()
  spectrum[["vects"]] <- vects
  spectrum[["vals"]] <- vals

  return(spectrum)
}

#' Build a Laplacian matrix
#'
#' Build a Laplacian matrix (combinatorial Laplacian or random-walk Laplacian)
#' from a symmetric (weighted) graph adjacency matrix.
#' @param adj_mat Symmetric adjacency matrix from which to build the Laplacian.
#' @param type_lap Type of Laplacian to build.
#' One of \code{"comb"} (combinatorial) or \code{"rw"} (random-walk).
#' @return The specified Laplacian matrix.
#' @examples
#' adj_mat <- matrix(c(1:9), nrow = 3)
#' build_laplacian(adj_mat, "rw")
#' @export

build_laplacian <- function(adj_mat, type_lap = c("comb", "rw")) {

  # check args
  type_lap <- match.arg(type_lap)
  adj_mat <- drop0(adj_mat)

  # initialize parameters
  degs_adj_mat <- apply(adj_mat, 1, sum)
  n <- nrow(adj_mat)

  # combinatorial Laplacian
  if (type_lap == "comb") {
    degs_matrix <- diag(degs_adj_mat)
    L <-  degs_matrix - adj_mat
  }

  # random-walk Laplacian
  else if (type_lap == "rw") {

    if (!all(degs_adj_mat != 0)) {
      stop("row sums of adj_mat must be non-zero")
    }

    inv_degs_matrix <- diag(degs_adj_mat ^ (-1))
    L <-  diag(n) - inv_degs_matrix %*% adj_mat
  }

  L <- drop0(L)

  return(L)
}

#' Run Laplace embedding
#'
#' Run Laplace embedding on a symmetric (weighted) adjacency matrix
#' with a specified number of eigenvalues and eigenvectors.
#' @param adj_mat Symmetric adjacency matrix to be embedded.
#' @param num_eigs Number of eigenvalues and eigenvectors for the embedding.
#' @param type_lap Type of Laplacian for the embedding.
#' One of \code{"comb"} (combinatorial) or \code{"rw"} (random-walk).
#' @return A list with two entries:
#' \code{vals} contains the length-\code{num_eigs} vector
#' of the first few eigenvalues of the Laplacian,
#' and \code{vects} contains an \code{nrow(adj_mat)} by \code{num_eigs} matrix
#' of the associated eigenvectors.
#' @examples
#' adj_mat <- matrix(c(1:9), nrow = 3)
#' run_laplace_embedding(adj_mat, 2, "rw")
#' @export

run_laplace_embedding <- function(adj_mat, num_eigs,
                                  type_lap = c("comb", "rw")) {

  # check args
  stopifnot(all.equal(num_eigs, as.integer(num_eigs)))
  stopifnot(num_eigs > 0)
  type_lap <- match.arg(type_lap)

  # build and embed Laplacian
  laplacian <- build_laplacian(adj_mat, type_lap)
  spectrum <- get_first_eigs(laplacian, num_eigs)

  return(spectrum)
}

#' Run motif embedding
#'
#' Calculate a motif adjacency matrix for a given motif and motif type,
#' restrict it to its largest connected component,
#' and then run Laplace embedding with specified Laplacian type and
#' number of eigenvalues and eigenvectors.
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
#' @return A list with 7 entries:
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
#' }
#' @examples
#' adj_mat <- matrix(c(1:9), nrow = 3)
#' run_motif_embedding(adj_mat, "M1", "func")
#' @export

run_motif_embedding <- function(adj_mat, motif_name,
                       motif_type = c("struc", "func"),
                       mam_weight_type = c("unweighted", "mean", "product"),
                       mam_method = c("sparse", "dense"),
                       num_eigs = 2, type_lap = c("comb", "rw"),
                       restrict = TRUE) {

  # check args
  adj_mat <- drop0(adj_mat)
  stopifnot(motif_name %in% get_motif_names())
  motif_type <- match.arg(motif_type)
  stopifnot(all.equal(num_eigs, as.integer(num_eigs)))
  mam_weight_type <- match.arg(mam_weight_type)
  mam_method <- match.arg(mam_method)
  stopifnot(num_eigs > 0)
  type_lap <- match.arg(type_lap)
  stopifnot(typeof(restrict) == typeof(TRUE))

  # build motif adjacency matrix
  motif_adj_mat <- build_motif_adjacency_matrix(adj_mat, motif_name,
                     motif_type, mam_weight_type, mam_method)

  if (restrict) {

    # restrict to largest connected component
    comps <- get_largest_component(motif_adj_mat)
    adj_mat_comps <- adj_mat[comps, comps, drop = FALSE]
    motif_adj_mat_comps <- motif_adj_mat[comps, comps, drop = FALSE]

    # Laplace embedding restricted
    spect <- run_laplace_embedding(motif_adj_mat_comps, num_eigs, type_lap)
  }

  else {
    comps <- NULL
    adj_mat_comps <- NULL
    motif_adj_mat_comps <- NULL

    # Laplace embedding unrestricted
    spect <- run_laplace_embedding(motif_adj_mat, num_eigs, type_lap)
  }

  # return list
  embedding <- list()
  embedding$adj_mat <- adj_mat
  embedding$motif_adj_mat <- motif_adj_mat
  embedding$comps <- comps
  embedding$adj_mat_comps <- adj_mat_comps
  embedding$motif_adj_mat_comps <- motif_adj_mat_comps
  embedding$vals <- spect$vals
  embedding$vects <- spect$vects

  return(embedding)
}
