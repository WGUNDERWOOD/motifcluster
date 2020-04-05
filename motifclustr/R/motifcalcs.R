# TODO remove run_motif_adj_calcs

#' Run motif adjacency matrix calculations
#'
#' Perform the matrix operations required to build a motif
#' adjacency matrix from adjacency and indicator matrices.
#' Parameter names default to building functional motif adjacency matrices.
#' @param G The original adjacency matrix.
#' Replaced by the single-edge adjacency matrix Gs for structural motifs.
#' @param J The directed indicator matrix.
#' Replaced by the single-edge indicator matrix Js for structural motifs.
#' @param Jn is the vertex-distinct indicator matrix.
#' Replaced by the missing-edge indicator matrix J0 for structural motifs.
#' @param Jd The double-edge indicator matrix.
#' @param Gd The double-edge adjacency matrix.
#' @param motif_name The motif to use for the motif adjacency matrix.
#' @return A motif adjacency matrix.
#' @keywords motif adjacency matrix calculation operation

run_motif_adjacency_calcs <- function(ind_mats, motif_name, motif_type, weight_type, method){

  if(weight_type == "unweighted"){
    # product calcs
  }
  if(weight_type == "product"){
    # product calcs
  }
  if(weight_type == "mean"){
    # mean calcs
  }






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
# TODO docs

mam_M1 <- function(adj_mat, motif_type, weight_type){

  if(weight_type == "unweighted"){
    if(motif_type == "func"){
      J = build_J(adj_mat)
      C <- t(J)*(J%*%J)
      return(C + t(C))
    }

    if(motif_type == "struc"){
      Js = build_Js(adj_mat)
      C <- t(Js)*(Js%*%Js)
      return(C + t(C))
    }
  }

  if(weight_type == "mean"){
    if(motif_type == "func"){
      J = build_J(adj_mat)
      G = adj_mat
      C <- t(J)*(J%*%G) + t(J)*(G%*%J) + t(G)*(J%*%J)
      return((C + t(C)) / 3)
    }

    if(motif_type == "struc"){
      Js = build_Js(adj_mat)
      Gs = build_Gs(adj_mat)
      C <- t(Js)*(Js%*%Gs) + t(Js)*(Gs%*%Js) + t(Gs)*(Js%*%Js)
      return((C + t(C)) / 3)
    }
  }

  if(weight_type == "product"){
    if(motif_type == "func"){
      G = adj_mat
      C <- t(G)*(G%*%G)
      return(C + t(C))
    }

    if(motif_type == "struc"){
      Gs = build_Gs(adj_mat)
      C <- t(Gs)*(Gs%*%Gs)
      return(C + t(C))
    }
  }
}

mam_Ms <- function(adj_mat, motif_type, weight_type){

  if(weight_type == "unweighted"){
    if(motif_type == "func"){
      J = build_J(adj_mat)
      return(J + t(J))
    }

    if(motif_type == "struc"){
      Js = build_Js(adj_mat)
      return(Js + t(Js))
    }
  }

  if(weight_type == "mean"){
    if(motif_type == "func"){
      G = adj_mat
      return(G + t(G))
    }

    if(motif_type == "struc"){
      Gs = build_Gs(adj_mat)
      return(Gs + t(Gs))
    }
  }

  if(weight_type == "product"){
    if(motif_type == "func"){
      G = adj_mat
      return(G + t(G))
    }

    if(motif_type == "struc"){
      Gs = build_Gs(adj_mat)
      return(Gs + t(Gs))
    }
  }
}

mam_Md <- function(adj_mat, weight_type){

  if(weight_type == "unweighted"){
    Jd = build_Jd(adj_mat)
    return(Jd)
  }

  if(weight_type == "mean"){
    Gd = build_Gd(adj_mat)
    return(Gd/2)
  }

  if(weight_type == "product"){
    G = adj_mat
    return(G * t(G))
  }
}
