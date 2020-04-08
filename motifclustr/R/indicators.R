#' Build sparse adjacency matrix
#'
#' Build the sparse adjacency matrix G from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return The adjacency matrix G in sparse form.

build_G <- function(adj_mat){
  G = drop0(adj_mat)
  return(G)
}

#' Build directed indicator matrix
#'
#' Build the sparse directed indicator matrix J from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A directed indicator matrix J in sparse form.

build_J <- function(adj_mat){
  J = drop0_killdiag(1*(adj_mat > 0))
  return(J)
}

#' Build single-edge adjacency matrix
#'
#' Build the sparse single-edge adjacency matrix Gs from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A single-edge adjacency matrix Gs in sparse form.

build_Gs <- function(adj_mat){
  J = build_J(adj_mat)
  Gs = drop0_killdiag(adj_mat*(1 - t(J)))
  return(Gs)
}

#' Build single-edge indicator matrix
#'
#' Build the sparse single-edge indicator matrix Js from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A single-edge indicator matrix Js in sparse form.

build_Js <- function(adj_mat){
  Gs = build_Gs(adj_mat)
  Js = drop0_killdiag(1*(Gs > 0))
  return(Js)
}

#' Build double-edge adjacency matrix
#'
#' Build the sparse double-edge adjacency matrix Gd from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A double-edge adjacency matrix Gd in sparse form.

build_Gd <- function(adj_mat){
  J = build_J(adj_mat)
  Gd = drop0_killdiag((adj_mat+t(adj_mat)) * J * t(J))
  return(Gd)
}

#' Build double-edge indicator matrix
#'
#' Build the sparse double-edge indicator matrix Jd from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A double-edge indicator matrix Jd in sparse form.

build_Jd <- function(adj_mat){
  Gd = build_Gd(adj_mat)
  Jd = drop0_killdiag(1*(Gd > 0))
  return(Jd)
}

#' Build missing-edge indicator matrix
#'
#' Build the missing-edge indicator matrix J0 from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A missing-edge indicator matrix J0.

build_J0 <- function(adj_mat){
  J0 = drop0_killdiag(1*((adj_mat+t(adj_mat)) == 0))
  return(J0)
}

#' Build vertex-distinct indicator matrix
#'
#' Build the vertex-distinct indicator matrix Jn from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A vertex-distinct indicator matrix Jn.

build_Jn <- function(adj_mat){
  Jn = drop0_killdiag(1+0*adj_mat)
  return(Jn)
}

#' Build identity matrix
#'
#' Build the sparse identity matrix Id from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return An identity matrix Id in sparse form.
#' @importFrom Matrix Diagonal

build_Id <- function(adj_mat){
  Id <- Diagonal(nrow(adj_mat))
  return(Id)
}

#' Build edge-and-diagonal indicator matrix
#'
#' Build the sparse edge-and-diagonal indicator matrix Je from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return An edge-and-diagonal indicator matrix Je in sparse form.

build_Je <- function(adj_mat){
  Id = build_Id(adj_mat)
  Je = drop0(1*(Id + ((adj_mat+t(adj_mat)) > 0)))
  return(Je)
}

#' Build product adjacency matrix
#'
#' Build the sparse product adjacency matrix Jp from a graph adjacency matrix.
#' @param adj_mat Original adjacency matrix.
#' @return A product adjacency matrix Jp in sparse form.

build_Gp <- function(adj_mat){
  Gp = drop0_killdiag(adj_mat*t(adj_mat))
  return(Gp)
}
