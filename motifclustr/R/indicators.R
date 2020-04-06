# TODO redo ind mat docs

#' Build indicator matrices
#'
#' Build the indicator matrices and variants of
#' adjacency matrices from a graph adjacency matrix.
#' @param adj_mat Adjacency matrix from which to build the indicator matrices.
#' @return A list with 8 entries:
#' G is the original adjacency matrix;
#' J is the directed indicator matrix;
#' J0 is the missing-edge indicator matrix;
#' Jn is the vertex-distinct indicator matrix;
#' Gs is the single-edge adjacency matrix;
#' Gd is the double-edge adjacency matrix
#' Js is the single-edge indicator matrix;
#' Jd is the double-edge indicator matrix.
#' @keywords indicator adjacency matrix

build_J <- function(adj_mat){
  J = drop0_killdiag(1*(adj_mat > 0))
  return(J)
}

build_Gs <- function(adj_mat){
  J = build_J(adj_mat)
  Gs = drop0_killdiag(adj_mat*(1 - t(J)))
  return(Gs)
}

build_Gd <- function(adj_mat){
  J = build_J(adj_mat)
  Gd = drop0_killdiag((adj_mat+t(adj_mat)) * J * t(J))
  return(Gd)
}

build_Js <- function(adj_mat){
  Gs = build_Gs(adj_mat)
  Js = drop0_killdiag(1*(Gs > 0))
  return(Js)
}

build_Jd <- function(adj_mat){
  Gd = build_Gd(adj_mat)
  Jd = drop0_killdiag(1*(Gd > 0))
  return(Jd)
}

build_J0 <- function(adj_mat){
  J0 = drop0_killdiag(1*((adj_mat+t(adj_mat)) == 0))
  return(J0)
}

build_Jn <- function(adj_mat){
  Jn = drop0_killdiag(1+0*adj_mat)
  return(Jn)
}

build_Id <- function(adj_mat){
  Id <- Diagonal(nrow(adj_mat))
  return(Id)
}

build_Je <- function(adj_mat){
  Je = drop0(1*((adj_mat+t(adj_mat)) > 0))
}
