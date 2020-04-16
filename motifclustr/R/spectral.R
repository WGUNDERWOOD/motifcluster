#' Compute first few eigenvalues and eigenvectors of a matrix
#'
#' Compute the first few eigenvalues (by magnitude) and
#' associated eigenvectors of a matrix.
#' @param mat Symmetric matrix for which eigenvalues and eigenvectors
#' are to be calculated.
#' @param num_eigs Number of eigenvalues and eigenvectors to calculate.
#' @return A list with two entries:
#' \code{vals} contains a length-\code{num_eigs} vector of the first few
#' eigenvalues,
#' and vects contains an \code{nrow(mat)} by \code{num_eigs} matrix
#' of the associated eigenvectors.
#' @importFrom RSpectra eigs
#' @keywords internal

get_first_eigs <- function(mat, num_eigs) {

  # TODO check returned data types.

  # check args
  if (!all.equal(num_eigs, as.integer(num_eigs))) {
    stop("num_eigs must be an integer.")
  }
  if (!(num_eigs > 0)) {
    stop("num_eigs must be at least 1.")
  }

  # get spectrum
  ans_eigs <- eigs(mat, num_eigs, which = "SM")

  # order eigenvalues and eigenvectors
  inds <- seq(num_eigs, 1, -1)
  vects <- Re(ans_eigs[["vectors"]])[, inds]
  vals <- Re(ans_eigs[["values"]])[inds]

  # return a list
  ans_spect <- list()
  ans_spect[["vects"]] <- vects
  ans_spect[["vals"]] <- vals

  return(ans_spect)
}

#' Build a Laplacian matrix
#'
#' Build a Laplacian matrix (combinatorial Laplacian or random-walk Laplacian)
#' from a symmetric (weighted) graph adjacency matrix.
#' @param adj_mat Symmetric adjacency matrix from which to build the Laplacian.
#' @param type_lap Type of Laplacian to build.
#' One of \code{"comb"} (combinatorial) or \code{"rw"} (random-walk).
#' @return The specified Laplacian matrix.
#' @export
# TODO examples

build_laplacian <- function(adj_mat, type_lap = c("comb", "rw")) {

  # check args
  type_lap <- match.arg(type_lap)

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
#' @export
# TODO examples

run_laplace_embedding <- function(adj_mat, num_eigs,
                                  type_lap = c("comb", "rw")) {

  # check args
  if (!all.equal(num_eigs, as.integer(num_eigs))) {
    stop("num_eigs must be an integer.")
  }
  if (!(num_eigs > 0)) {
    stop("num_eigs must be at least 1.")
  }
  type_lap <- match.arg(type_lap)

  # build and embed Laplacian
  laplacian <- build_laplacian(adj_mat, type_lap)
  ans_spect <- get_first_eigs(laplacian, num_eigs)

  return(ans_spect)
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
#' @param type_lap Type of Laplacian for the embedding.
#' One of \code{"comb"} or \code{"rw"}.
#' @param num_eigs Number of eigenvalues and eigenvectors for the embedding.
#' @return A list with 7 entries:
#' \itemize{
#'   \item \code{adj_mat}: the original adjacency matrix.
#'   \item \code{motif_adj_mat}: the motif adjacency matrix.
#'   \item \code{comps}: the indices of the largest connected component
#'     of the motif adjacency matrix.
#'   \item \code{adj_mat_comps}: the original adjacency matrix restricted
#'     to the largest connected component of the motif adjacency matrix.
#'   \item \code{motif_adj_mat_comps}: the motif adjacency matrix restricted
#'     to its largest connected component.
#'   \item \code{vals}: the eigenvalues associated with the Laplace embedding
#'     of the motif adjacency matrix.
#'   \item \code{vects}: the eigenvectors associated with the Laplace embedding
#'     of the motif adjacency matrix.
#' }
#' @export
# TODO examples

# TODO data types in vals and vects

run_motif_embedding <- function(adj_mat, motif_name,
                       motif_type = c("func", "struc"), num_eigs, type_lap) {

  # check args
  if (!(motif_name %in% get_motif_names())) {
    stop("Invalid motif name.")
  }
  motif_type <- match.arg(motif_type)
  if (!all.equal(num_eigs, as.integer(num_eigs))) {
    stop("num_eigs must be an integer.")
  }
  if (!(num_eigs > 0)) {
    stop("num_eigs must be at least 1.")
  }
  type_lap <- match.arg(type_lap)

  # build motif adjacency matrix
  motif_adj_mat <- build_motif_adjacency_matrix(adj_mat, motif_name,
                                                motif_type)

  # restrict to largest connected component
  comps <- get_largest_component(motif_adj_mat)
  adj_mat_comps <- adj_mat[comps, comps, drop = FALSE]
  motif_adj_mat_comps <- motif_adj_mat[comps, comps, drop = FALSE]

  # Laplace embedding
  spect <- run_laplace_embedding(motif_adj_mat, num_eigs, type_lap)

  # return list
  ans <- list()
  ans$adj_mat <- adj_mat
  ans$motif_adj_mat <- motif_adj_mat
  ans$comps <- comps
  ans$adj_mat_comps <- adj_mat_comps
  ans$motif_adj_mat_comps <- motif_adj_mat_comps
  ans$vals <- spect$vals
  ans$vects <- spect$vects

  return(ans)
}
