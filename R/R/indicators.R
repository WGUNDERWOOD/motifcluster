#' Build sparse adjacency matrix
#'
#' Build the sparse adjacency matrix \code{G} from a graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return The adjacency matrix \code{G} in sparse form.
#' @keywords internal

build_G <- function(adj_mat) {
  G <- drop0(adj_mat)
  return(G)
}

#' Build directed indicator matrix
#'
#' Build the sparse directed indicator matrix \code{J}
#' from a graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A directed indicator matrix \code{J} in sparse form.
#' @keywords internal

build_J <- function(adj_mat) {
  J <- drop0_killdiag(1 * (adj_mat > 0))
  return(J)
}

#' Build single-edge adjacency matrix
#'
#' Build the sparse single-edge adjacency matrix \code{Gs} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A single-edge adjacency matrix \code{Gs} in sparse form.
#' @keywords internal

build_Gs <- function(adj_mat) {
  J <- build_J(adj_mat)
  Gs <- drop0_killdiag(adj_mat - adj_mat * t(J))
  return(Gs)
}

#' Build single-edge indicator matrix
#'
#' Build the sparse single-edge indicator matrix \code{Js} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A single-edge indicator matrix \code{Js} in sparse form.
#' @keywords internal

build_Js <- function(adj_mat) {
  Gs <- build_Gs(adj_mat)
  Js <- drop0_killdiag(1 * (Gs > 0))
  return(Js)
}

#' Build double-edge adjacency matrix
#'
#' Build the sparse double-edge adjacency matrix \code{Gd} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A double-edge adjacency matrix \code{Gd} in sparse form.
#' @keywords internal

build_Gd <- function(adj_mat) {
  J <- build_J(adj_mat)
  Gd <- drop0_killdiag((adj_mat + t(adj_mat)) * J * t(J))
  return(Gd)
}

#' Build double-edge indicator matrix
#'
#' Build the sparse double-edge indicator matrix \code{Jd} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A double-edge indicator matrix \code{Jd} in sparse form.
#' @keywords internal

build_Jd <- function(adj_mat) {
  Gd <- build_Gd(adj_mat)
  Jd <- drop0_killdiag(1 * (Gd > 0))
  return(Jd)
}

#' Build missing-edge indicator matrix
#'
#' Build the missing-edge indicator matrix \code{J0} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A missing-edge indicator matrix \code{J0}.
#' @keywords internal

build_J0 <- function(adj_mat) {
  J0 <- drop0_killdiag(1 * ((adj_mat + t(adj_mat)) == 0))
  return(J0)
}

#' Build vertex-distinct indicator matrix
#'
#' Build the vertex-distinct indicator matrix \code{Jn} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A vertex-distinct indicator matrix \code{Jn}.
#' @keywords internal

build_Jn <- function(adj_mat) {
  Jn <- drop0_killdiag(1 + 0 * adj_mat)
  return(Jn)
}

#' Build identity matrix
#'
#' Build the sparse identity matrix \code{Id} from a graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return An identity matrix \code{Id} in sparse form.
#' @keywords internal
#' @importFrom Matrix Diagonal

build_Id <- function(adj_mat) {
  Id <- Diagonal(nrow(adj_mat))
  return(Id)
}

#' Build edge-and-diagonal indicator matrix
#'
#' Build the sparse edge-and-diagonal indicator matrix \code{Je} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return An edge-and-diagonal indicator matrix \code{Je} in sparse form.
#' @keywords internal

build_Je <- function(adj_mat) {
  Id <- build_Id(adj_mat)
  Je <- drop0(1 * (Id + ((adj_mat + t(adj_mat)) > 0)))
  return(Je)
}

#' Build product adjacency matrix
#'
#' Build the sparse product adjacency matrix \code{Jp} from a
#' graph adjacency matrix.
#' @param adj_mat The original adjacency matrix.
#' @return A product adjacency matrix \code{Jp} in sparse form.
#' @keywords internal

build_Gp <- function(adj_mat) {
  Gp <- drop0_killdiag(adj_mat * t(adj_mat))
  return(Gp)
}
