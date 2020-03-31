#' Build a motif adjacency matrix
#'
#' Build a motif adjacency matrix from an adjacency matrix.
#' @param adj_mat Adjacency matrix from which to build the motif adjacency matrix.
#' @param motif_name Motif used for the motif adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @return A motif adjacency matrix.
#' @keywords motif adjacency matrix
#' @export
#' @examples
# TODO example

build_motif_adjacency_matrix <- function(adj_mat, motif_name, motif_type=c("func","struc")){

  # check args
  if(!is.matrix(adj_mat)){
    stop("adj_mat must be a matrix.")
  }
  if(!(motif_name %in% get_motif_names())){
    stop("Invalid motif name.")
  }
  motif_type <- match.arg(motif_type)

  # initialize parameters
  G <- adj_mat
  IM <- build_indicator_matrices(G)

  # build functional motif adjacency matrix
  if(motif_type=="func"){
    # TODO redo the args below
    #motif_adj_mat <- run_motif_adjacency_calcs(G, G, IM$Gd, IM$J, IM$Jn, IM$J, IM$Jd, motif_name)
  }

  # build structural motif adjacency matrix
  else if(motif_type=="struc"){
    # TODO redo the args below
    #motif_adj_mat <- run_motif_adjacency_calcs(G, IM$Gs, IM$Gd, IM$J, IM$J0, IM$Js, IM$Jd, motif_name)
  }

  motifadj <- unname(drop0(motifadj))

  return(motifadj)
}

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
#' @export
#' @examples
# TODO example

build_indicator_matrices <- function(adj_mat){

  G  <- adj_mat
  J  <- drop0_killdiag( 1*(G > 0) )
  J0 <- drop0_killdiag( 1*((G+t(G)) == 0 ) )
  Jn <- drop0_killdiag( 1+0*G )
  Gs <- drop0_killdiag( G*(1 - t(J)) )
  Gd <- drop0_killdiag( (G+t(G)) * J * t(J) )
  Js <- drop0_killdiag( 1*(Gs > 0) )
  Jd <- drop0_killdiag( 1*(Gd > 0) )

  return(list(G=G, J=J, J0=J0, Jn=Jn, Gs=Gs, Gd=Gd, Js=Js, Jd=Jd))
}

#' Run motif adjacency matrix calculations
#'
#' Perform the matrix operations required to build a motif
#' adjacency matrix from adjacency and indicator matrices.
#' @param G The original adjacency matrix.
(re)
#' @param J The directed indicator matrix.
#' @param Jn is the vertex-distinct indicator matrix.
#' @param Jd The double-edge indicator matrix.
#' @param Gd The double-edge adjacency matri.
# TODO finish these docs

run_motif_adjacency_calcs <- function(G, J, Jn, Jd, Gd, motif_name){

  if(motif_name == "Ms"){
    motif_adj_mat <- G + t(G)
  }

  else if(motif_name == "Md"){
    motif_adj_mat <- Gd / 2
  }

  else if(motif_name == "M1"){
    C <- t(J)*(J%*%G) + t(J)*(G%*%J) + t(G)*(J%*%J)
    motif_adj_mat <- (C + t(C)) / 3
  }

  else if(motif_name == "M2"){
    C <- t(J)*(Jd%*%G) + t(J)*(Gd%*%J) + t(G)*(Jd%*%J)
    C <- C + t(J)*(J%*%Gd) + t(J)*(G%*%Jd) + t(G)*(J%*%Jd)
    C <- C + Jd*(J%*%G) + Jd*(G%*%J) + Gd*(J%*%J)
    motif_adj_mat <- (C + t(C)) / 4
  }

  else if(motif_name == "M3"){
    C <- J*(Jd%*%Gd) + J*(Gd%*%Jd) + G*(Jd%*%Jd)
    C <- C + Jd*(Jd%*%G) + Jd*(Gd%*%J) + Gd*(Jd%*%J)
    C <- C + Jd*(J%*%Gd) + Jd*(G%*%Jd) + Gd*(J%*%Jd)
    motif_adj_mat <- (C + t(C)) / 5
  }

  else if(motif_name == "M4"){
    motif_adj_mat <- (Jd*(Jd%*%Gd) + Jd*(Gd%*%Jd) + Gd*(Jd%*%Jd)) / 6
  }

  else if(motif_name == "M5"){
    C <- J*(J%*%G) + J*(G%*%J) + G*(J%*%J)
    C <- C + J*(J%*%t(G)) + J*(G%*%t(J)) + G*(J%*%t(J))
    C <- C + J*(t(J)%*%G) + J*(t(G)%*%J) + G*(t(J)%*%J)
    motif_adj_mat <- (C + t(C)) / 3
  }

  else if(motif_name == "M6"){
    C <- J*(J%*%Gd) + J*(G%*%Jd) + G*(J%*%Jd)
    Cprime <- Jd*(t(J)%*%G) + Jd*(t(G)%*%J) + Gd*(t(J)%*%J)
    motif_adj_mat <- (C + t(C) + Cprime) / 4
  }

  else if(motif_name == "M7"){
    C <- J*(Jd%*%G) + J*(Gd%*%J) + G*(Jd%*%J)
    Cprime <- Jd*(J%*%t(G)) + Jd*(G%*%t(J)) + Gd*(J%*%t(J))
    motif_adj_mat <- (C + t(C) + Cprime) / 4
  }

  else if(motif_name == "M8"){
    C <- J*(G%*%Jn) + G*(J%*%Jn)
    Cprime <- Jn*(t(J)%*%G) + Jn*(t(G)%*%J)
    motif_adj_mat <- (C + t(C) + Cprime) / 2
  }

  else if(motif_name == "M9"){
    C <- J*(Jn%*%t(G)) + G*(Jn%*%t(J))
    C <- C + Jn*(J%*%G) + Jn*(G%*%J)
    C <- C + J*(t(G)%*%Jn) + G*(t(J)%*%Jn)
    motif_adj_mat <- (C + t(C)) / 2
  }

  else if(motif_name == "M10"){
    C <- J*(Jn%*%G) + G*(Jn%*%J)
    Cprime <- Jn*(J%*%t(G)) + Jn*(G%*%t(J))
    motif_adj_mat <- (C + t(C) + Cprime) / 2
  }

  else if(motif_name == "M11"){
    C <- Jd*(G%*%Jn) + Gd*(J%*%Jn)
    C <- C + Jn*(Jd%*%G) + Jn*(Gd%*%J)
    C <- C + J*(Gd%*%Jn) + G*(Jd%*%Jn)
    motif_adj_mat <- (C + t(C)) / 3
  }

  else if(motif_name == "M12"){
    C <- Jd*(Jn%*%G) + Gd*(Jn%*%J)
    C <- C + Jn*(J%*%Gd) + Jn*(G%*%Jd)
    C <- C + J*(Jn%*%Gd) + G*(Jn%*%Jd)
    motif_adj_mat <- (C + t(C)) / 3
  }

  else if(motif_name == "M13"){
    C <- Jd*(Gd%*%Jn) + Gd*(Jd%*%Jn) + Jn*(Jd%*%Gd)
    motif_adj_mat <- (C + t(C)) / 4
  }

  else if(motif_name == "coll"){
    C <- Jn*(J%*%t(G)) + Jn*(G%*%t(J))
    motif_adj_mat <- C / 2
  }

  else if(motif_name == "expa"){
    C <- Jn*(t(J)%*%G) + Jn*(t(G)%*%J)
    motif_adj_mat <- C / 2
  }

  else if(motif_name == "path"){
    C <- Jn*(J%*%G) + Jn*(G%*%J)
    motif_adj_mat <- (C + t(C)) / 2
  }

  return(motif_adj_mat)
}

drop0_killdiag <- function(mat){

  # Set the diagonal entries of a matrix to zero,
  # and convert it to sparse form.

  ans <- mat
  diag(ans) <- 0
  ans <- drop0(ans)

  return(ans)
}

largestComponent <- function(G){

  # Return the vertices in the largest connected component
  # associated with an adjacency matrix.

  n <- nrow(G)
  Gr <- graph_from_adjacency_matrix(G, weighted=TRUE)
  comps <- components(Gr)
  verts_to_keep <- (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)
}

get_motif_names <- function()
