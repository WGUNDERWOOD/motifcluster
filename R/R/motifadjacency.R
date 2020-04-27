#' Build a motif adjacency matrix
#'
#' Build a motif adjacency matrix from an adjacency matrix.
#'
#' Entry (\emph{i}, \emph{j}) of a motif adjacency matrix is the
#' sum of the weights of all motifs containing both
#' nodes \emph{i} and \emph{j}.
#' The motif is specified by name and the type of motif instance can be one of:
#' \itemize{
#'   \item Functional: motifs should appear as subgraphs.
#'   \item Structural: motifs should appear as induced subgraphs.
#' }
#' The weighting scheme can be one of:
#' \itemize{
#'   \item Unweighted: the weight of any motif instance is one.
#'   \item Mean: the weight of any motif instance
#'     is the mean of its edge weights.
#'   \item Product: the weight of any motif instance
#'     is the product of its edge weights.
#' }
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_name Motif used for the motif adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' One of \code{"func"} or \code{"struc"}.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' The sparse formulation avoids generating large dense matrices
#' so tends to be faster for large sparse graphs.
#' @return A motif adjacency matrix.
#' @importFrom Matrix drop0 t
#' @examples
#' adj_mat <- matrix(c(1:9), nrow = 3)
#' build_motif_adjacency_matrix(adj_mat, "M1", "func", "mean")
#' @export

build_motif_adjacency_matrix <- function(adj_mat, motif_name,
  motif_type = c("struc", "func"),
  mam_weight_type = c("unweighted", "mean", "poisson"),
  mam_method = c("sparse", "dense")) {

  # check args
  stopifnot(motif_name %in% get_motif_names())
  motif_type <- match.arg(motif_type)
  mam_weight_type <- match.arg(mam_weight_type, mam_weight_type)
  mam_method <- match.arg(mam_method)

  if (motif_name == "Ms") {
    return(mam_Ms(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "Md") {
    return(mam_Md(adj_mat, mam_weight_type))
  }

  if (motif_name == "M1") {
    return(mam_M1(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M2") {
    return(mam_M2(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M3") {
    return(mam_M3(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M4") {
    return(mam_M4(adj_mat, mam_weight_type))
  }

  if (motif_name == "M5") {
    return(mam_M5(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M6") {
    return(mam_M6(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M7") {
    return(mam_M7(adj_mat, motif_type, mam_weight_type))
  }

  if (motif_name == "M8") {
    return(mam_M8(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "M9") {
    return(mam_M9(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "M10") {
    return(mam_M10(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "M11") {
    return(mam_M11(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "M12") {
    return(mam_M12(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "M13") {
    return(mam_M13(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "Mcoll") {
    return(mam_Mcoll(adj_mat, motif_type, mam_weight_type, mam_method))
  }

  if (motif_name == "Mexpa") {
    return(mam_Mexpa(adj_mat, motif_type, mam_weight_type, mam_method))
  }
}

#' Perform the motif adjacency matrix calculations for motif Ms
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_Ms <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      return(J + t(J))
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      return(Js + t(Js))
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      return(G + t(G))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      return(Gs + t(Gs))
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      return(G + t(G))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      return(Gs + t(Gs))
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif Md
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_Md <- function(adj_mat, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    Jd <- build_Jd(adj_mat)
    return(Jd)
  }

  if (mam_weight_type == "mean") {
    Gd <- build_Gd(adj_mat)
    return(Gd / 2)
  }

  if (mam_weight_type == "product") {
    Gp <- build_Gp(adj_mat)
    return(Gp)
  }
}

#' Perform the motif adjacency matrix calculations for motif M1
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M1 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      C <- t(J) * (J %*% J)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      C <- t(Js) * (Js %*% Js)
      return(C + t(C))
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      G <- build_G(adj_mat)
      C <- t(J) * (J %*% G) + t(J) * (G %*% J) + t(G) * (J %*% J)
      return((C + t(C)) / 3)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Gs <- build_Gs(adj_mat)
      C <- t(Js) * (Js %*% Gs) + t(Js) * (Gs %*% Js) + t(Gs) * (Js %*% Js)
      return((C + t(C)) / 3)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      C <- t(G) * (G %*% G)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      C <- t(Gs) * (Gs %*% Gs)
      return(C + t(C))
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M2
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M2 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- t(J) * (Jd %*% J) + t(J) * (J %*% Jd) + Jd * (J %*% J)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- t(Js) * (Jd %*% Js) + t(Js) * (Js %*% Jd) + Jd * (Js %*% Js)
      return(C + t(C))
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      G <- build_G(adj_mat)
      C <- t(J) * (Jd %*% G) + t(J) * (Gd %*% J) + t(G) * (Jd %*% J)
      C <- C + t(J) * (J %*% Gd) + t(J) * (G %*% Jd) + t(G) * (J %*% Jd)
      C <- C + Jd * (J %*% G) + Jd * (G %*% J) + Gd * (J %*% J)
      return((C + t(C)) / 4)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gs <- build_Gs(adj_mat)
      Gd <- build_Gd(adj_mat)
      C <- t(Js) * (Jd %*% Gs) + t(Js) * (Gd %*% Js) + t(Gs) * (Jd %*% Js)
      C <- C + t(Js) * (Js %*% Gd) + t(Js) * (Gs %*% Jd) + t(Gs) * (Js %*% Jd)
      C <- C + Jd * (Js %*% Gs) + Jd * (Gs %*% Js) + Gd * (Js %*% Js)
      return((C + t(C)) / 4)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- t(G) * (Gp %*% G) + t(G) * (G %*% Gp) + Gp * (G %*% G)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- t(Gs) * (Gp %*% Gs) + t(Gs) * (Gs %*% Gp) + Gp * (Gs %*% Gs)
      return(C + t(C))
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M3
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M3 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- J * (Jd %*% Jd) + Jd * (Jd %*% J) + Jd * (J %*% Jd)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- Js * (Jd %*% Jd) + Jd * (Jd %*% Js) + Jd * (Js %*% Jd)
      return(C + t(C))
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      G <- build_G(adj_mat)
      C <- J * (Jd %*% Gd) + J * (Gd %*% Jd) + G * (Jd %*% Jd)
      C <- C + Jd * (Jd %*% G) + Jd * (Gd %*% J) + Gd * (Jd %*% J)
      C <- C + Jd * (J %*% Gd) + Jd * (G %*% Jd) + Gd * (J %*% Jd)
      return((C + t(C)) / 5)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gs <- build_Gs(adj_mat)
      Gd <- build_Gd(adj_mat)
      C <- Js * (Jd %*% Gd) + Js * (Gd %*% Jd) + Gs * (Jd %*% Jd)
      C <- C + Jd * (Jd %*% Gs) + Jd * (Gd %*% Js) + Gd * (Jd %*% Js)
      C <- C + Jd * (Js %*% Gd) + Jd * (Gs %*% Jd) + Gd * (Js %*% Jd)
      return((C + t(C)) / 5)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- G * (Gp %*% Gp) + Gp * (Gp %*% G) + Gp * (G %*% Gp)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- Gs * (Gp %*% Gp) + Gp * (Gp %*% Gs) + Gp * (Gs %*% Gp)
      return(C + t(C))
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M4
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M4 <- function(adj_mat, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    Jd <- build_Jd(adj_mat)
    return(Jd * (Jd %*% Jd))
  }

  if (mam_weight_type == "mean") {
    Jd <- build_Jd(adj_mat)
    Gd <- build_Gd(adj_mat)
    return((Jd * (Jd %*% Gd) + Jd * (Gd %*% Jd) + Gd * (Jd %*% Jd)) / 6)
  }

  if (mam_weight_type == "product") {
    Gp <- build_Gp(adj_mat)
    return(Gp * (Gp %*% Gp))
  }
}

#' Perform the motif adjacency matrix calculations for motif M5
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M5 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      C <- J * (J %*% J) + J * (J %*% t(J)) + J * (t(J) %*% J)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      C <- Js * (Js %*% Js) + Js * (Js %*% t(Js)) + Js * (t(Js) %*% Js)
      return(C + t(C))
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      G <- build_G(adj_mat)
      C <- J * (J %*% G) + J * (G %*% J) + G * (J %*% J)
      C <- C + J * (J %*% t(G)) + J * (G %*% t(J)) + G * (J %*% t(J))
      C <- C + J * (t(J) %*% G) + J * (t(G) %*% J) + G * (t(J) %*% J)
      return((C + t(C)) / 3)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Gs <- build_Gs(adj_mat)
      C <- Js * (Js %*% Gs) + Js * (Gs %*% Js) + Gs * (Js %*% Js)
      C <- C + Js * (Js %*% t(Gs)) + Js * (Gs %*% t(Js)) + Gs * (Js %*% t(Js))
      C <- C + Js * (t(Js) %*% Gs) + Js * (t(Gs) %*% Js) + Gs * (t(Js) %*% Js)
      return((C + t(C)) / 3)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      C <- G * (G %*% G) + G * (G %*% t(G)) + G * (t(G) %*% G)
      return(C + t(C))
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      C <- Gs * (Gs %*% Gs) + Gs * (Gs %*% t(Gs)) + Gs * (t(Gs) %*% Gs)
      return(C + t(C))
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M6
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M6 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- J * (J %*% Jd)
      Cprime <- Jd * (t(J) %*% J)
      return(C + t(C) + Cprime)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- Js * (Js %*% Jd)
      Cprime <- Jd * (t(Js) %*% Js)
      return(C + t(C) + Cprime)
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      G <- build_G(adj_mat)
      C <- J * (J %*% Gd) + J * (G %*% Jd) + G * (J %*% Jd)
      Cprime <- Jd * (t(J) %*% G) + Jd * (t(G) %*% J) + Gd * (t(J) %*% J)
      return((C + t(C) + Cprime) / 4)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Gs <- build_Gs(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      C <- Js * (Js %*% Gd) + Js * (Gs %*% Jd) + Gs * (Js %*% Jd)
      Cprime <- Jd * (t(Js) %*% Gs) + Jd * (t(Gs) %*% Js) + Gd * (t(Js) %*% Js)
      return((C + t(C) + Cprime) / 4)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- G * (G %*% Gp)
      Cprime <- Gp * (t(G) %*% G)
      return(C + t(C) + Cprime)
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- Gs * (Gs %*% Gp)
      Cprime <- Gp * (t(Gs) %*% Gs)
      return(C + t(C) + Cprime)
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M7
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M7 <- function(adj_mat, motif_type, mam_weight_type) {

  if (mam_weight_type == "unweighted") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- J * (Jd %*% J)
      Cprime <- Jd * (J %*% t(J))
      return(C + t(C) + Cprime)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Jd <- build_Jd(adj_mat)
      C <- Js * (Jd %*% Js)
      Cprime <- Jd * (Js %*% t(Js))
      return(C + t(C) + Cprime)
    }
  }

  if (mam_weight_type == "mean") {
    if (motif_type == "func") {
      J <- build_J(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      G <- build_G(adj_mat)
      C <- J * (Jd %*% G) + J * (Gd %*% J) + G * (Jd %*% J)
      Cprime <- Jd * (J %*% t(G)) + Jd * (G %*% t(J)) + Gd * (J %*% t(J))
      return((C + t(C) + Cprime) / 4)
    }

    if (motif_type == "struc") {
      Js <- build_Js(adj_mat)
      Gs <- build_Gs(adj_mat)
      Jd <- build_Jd(adj_mat)
      Gd <- build_Gd(adj_mat)
      C <- Js * (Jd %*% Gs) + Js * (Gd %*% Js) + Gs * (Jd %*% Js)
      Cprime <- Jd * (Js %*% t(Gs)) + Jd * (Gs %*% t(Js)) + Gd * (Js %*% t(Js))
      return((C + t(C) + Cprime) / 4)
    }
  }

  if (mam_weight_type == "product") {
    if (motif_type == "func") {
      G <- build_G(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- G * (Gp %*% G)
      Cprime <- Gp * (G %*% t(G))
      return(C + t(C) + Cprime)
    }

    if (motif_type == "struc") {
      Gs <- build_Gs(adj_mat)
      Gp <- build_Gp(adj_mat)
      C <- Gs * (Gp %*% Gs)
      Cprime <- Gp * (Gs %*% t(Gs))
      return(C + t(C) + Cprime)
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M8
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M8 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- J * (J %*% Jn)
        Cprime <- Jn * (t(J) %*% J)
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (Js %*% J0)
        Cprime <- J0 * (t(Js) %*% Js)
        return(C + t(C) + Cprime)
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- J * (G %*% Jn) + G * (J %*% Jn)
        Cprime <- Jn * (t(J) %*% G) + Jn * (t(G) %*% J)
        return((C + t(C) + Cprime) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (Gs %*% J0) + Gs * (Js %*% J0)
        Cprime <- J0 * (t(Js) %*% Gs) + J0 * (t(Gs) %*% Js)
        return((C + t(C) + Cprime) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- G * (G %*% Jn)
        Cprime <- Jn * (t(G) %*% G)
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Gs * (Gs %*% J0)
        Cprime <- J0 * (t(Gs) %*% Gs)
        return(C + t(C) + Cprime)
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(J, J) - J * J
        Cprime <- t(J) %*% J - Id * (t(J) %*% J)
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Js, Js) - Js * (Js %*% Je)
        Cprime <- t(Js) %*% Js - Je * (t(Js) %*% Js)
        return(drop0(C + t(C) + Cprime))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(J, G) - J * G + a_b_one(G, J) - G * J
        Cprime <- t(J) %*% G - Id * (t(J) %*% G)
        Cprime <- Cprime + t(G) %*% J - Id * (t(G) %*% J)
        return(drop0(C + t(C) + Cprime) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Js, Gs) - Js * (Gs %*% Je)
        C <- C + a_b_one(Gs, Js) - Gs * (Js %*% Je)
        Cprime <- t(Js) %*% Gs - Je * (t(Js) %*% Gs)
        Cprime <- Cprime + t(Gs) %*% Js - Je * (t(Gs) %*% Js)
        return(drop0(C + t(C) + Cprime) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(G, G) - G * G
        Cprime <- t(G) %*% G - Id * (t(G) %*% G)
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Gs, Gs) - Gs * (Gs %*% Je)
        Cprime <- t(Gs) %*% Gs - Je * (t(Gs) %*% Gs)
        return(drop0(C + t(C) + Cprime))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M9
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M9 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- J * (Jn %*% t(J)) + Jn * (J %*% J) + J * (t(J) %*% Jn)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (J0 %*% t(Js)) + J0 * (Js %*% Js) + Js * (t(Js) %*% J0)
        return(C + t(C))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- J * (Jn %*% t(G)) + G * (Jn %*% t(J))
        C <- C + Jn * (J %*% G) + Jn * (G %*% J)
        C <- C + J * (t(G) %*% Jn) + G * (t(J) %*% Jn)
        return((C + t(C)) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (J0 %*% t(Gs)) + Gs * (J0 %*% t(Js))
        C <- C + J0 * (Js %*% Gs) + J0 * (Gs %*% Js)
        C <- C + Js * (t(Gs) %*% J0) + Gs * (t(Js) %*% J0)
        return((C + t(C)) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- G * (Jn %*% t(G)) + Jn * (G %*% G) + G * (t(G) %*% Jn)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Gs * (J0 %*% t(Gs)) + J0 * (Gs %*% Gs) + Gs * (t(Gs) %*% J0)
        return(C + t(C))
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(J, t(J)) - 2 * J * t(J) + J %*% J
        C <- C - Id * (J %*% J) + a_b_one(J, t(J))
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Js, t(Js)) - Js * (Je %*% t(Js))
        C <- C + Js %*% Js - Je * (Js %*% Js)
        C <- C + a_b_one(Js, t(Js)) - Js * (t(Js) %*% Je)
        return(drop0(C + t(C)))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(J, t(G)) - 2 * J * t(G) + J %*% G
        C <- C + a_one_b(G, t(J)) - 2 * G * t(J) + G %*% J
        C <- C - Id * (J %*% G) + a_b_one(J, t(G))
        C <- C - Id * (G %*% J) + a_b_one(G, t(J))
        return(drop0(C + t(C)) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Js, t(Gs)) - Js * (Je %*% t(Gs))
        C <- C + a_one_b(Gs, t(Js)) - Gs * (Je %*% t(Js))
        C <- C + Js %*% Gs - Je * (Js %*% Gs)
        C <- C + a_b_one(Js, t(Gs)) - Js * (t(Gs) %*% Je)
        C <- C + Gs %*% Js - Je * (Gs %*% Js)
        C <- C + a_b_one(Gs, t(Js)) - Gs * (t(Js) %*% Je)
        return(drop0(C + t(C)) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(G, t(G)) - 2 * G * t(G) + G %*% G
        C <- C - Id * (G %*% G) + a_b_one(G, t(G))
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Gs, t(Gs)) - Gs * (Je %*% t(Gs))
        C <- C + Gs %*% Gs - Je * (Gs %*% Gs)
        C <- C + a_b_one(Gs, t(Gs)) - Gs * (t(Gs) %*% Je)
        return(drop0(C + t(C)))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M10
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M10 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- J * (Jn %*% J)
        Cprime <- Jn * (J %*% t(J))
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (J0 %*% Js)
        Cprime <- J0 * (Js %*% t(Js))
        return(C + t(C) + Cprime)
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- J * (Jn %*% G) + G * (Jn %*% J)
        Cprime <- Jn * (J %*% t(G)) + Jn * (G %*% t(J))
        return((C + t(C) + Cprime) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Js * (J0 %*% Gs) + Gs * (J0 %*% Js)
        Cprime <- J0 * (Js %*% t(Gs)) + J0 * (Gs %*% t(Js))
        return((C + t(C) + Cprime) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- G * (Jn %*% G)
        Cprime <- Jn * (G %*% t(G))
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Gs * (J0 %*% Gs)
        Cprime <- J0 * (Gs %*% t(Gs))
        return(C + t(C) + Cprime)
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(J, J) - J * J
        Cprime <- J %*% t(J) - Id * (J %*% t(J))
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Js, Js) - Js * (Je %*% Js)
        Cprime <- Js %*% t(Js) - Je * (Js %*% t(Js))
        return(drop0(C + t(C) + Cprime))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(J, G) - J * G + a_one_b(G, J) - G * J
        Cprime <- J %*% t(G) - Id * (J %*% t(G))
        Cprime <- Cprime + G %*% t(J) - Id * (G %*% t(J))
        return(drop0(C + t(C) + Cprime) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Js, Gs) - Js * (Je %*% Gs)
        C <- C + a_one_b(Gs, Js) - Gs * (Je %*% Js)
        Cprime <- Js %*% t(Gs) - Je * (Js %*% t(Gs))
        Cprime <- Cprime + Gs %*% t(Js) - Je * (Gs %*% t(Js))
        return(drop0(C + t(C) + Cprime) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_one_b(G, G) - G * G
        Cprime <- G %*% t(G) - Id * (G %*% t(G))
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_one_b(Gs, Gs) - Gs * (Je %*% Gs)
        Cprime <- Gs %*% t(Gs) - Je * (Gs %*% t(Gs))
        return(drop0(C + t(C) + Cprime))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M11
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M11 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Jn <- build_Jn(adj_mat)
        J <- build_J(adj_mat)
        C <- Jd * (J %*% Jn) + Jn * (Jd %*% J) + J * (Jd %*% Jn)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        J0 <- build_J0(adj_mat)
        Js <- build_Js(adj_mat)
        C <- Jd * (Js %*% J0) + J0 * (Jd %*% Js) + Js * (Jd %*% J0)
        return(C + t(C))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Jn <- build_Jn(adj_mat)
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        C <- Jd * (G %*% Jn) + Gd * (J %*% Jn)
        C <- C + Jn * (Jd %*% G) + Jn * (Gd %*% J)
        C <- C + J * (Gd %*% Jn) + G * (Jd %*% Jn)
        return((C + t(C)) / 3)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        J0 <- build_J0(adj_mat)
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- Jd * (Gs %*% J0) + Gd * (Js %*% J0)
        C <- C + J0 * (Jd %*% Gs) + J0 * (Gd %*% Js)
        C <- C + Js * (Gd %*% J0) + Gs * (Jd %*% J0)
        return((C + t(C)) / 3)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- Gp * (G %*% Jn) + Jn * (Gp %*% G) + G * (Gp %*% Jn)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        J0 <- build_J0(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- Gp * (Gs %*% J0) + J0 * (Gp %*% Gs) + Gs * (Gp %*% J0)
        return(C + t(C))
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Id <- build_Id(adj_mat)
        J <- build_J(adj_mat)
        C <- a_b_one(Jd, J) - Jd * J
        C <- C + Jd %*% J - Id * (Jd %*% J)
        C <- C + a_b_one(J, Jd) - J * Jd
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Je <- build_Je(adj_mat)
        Js <- build_Js(adj_mat)
        C <- a_b_one(Jd, Js) - Jd * (Js %*% Je)
        C <- C + Jd %*% Js - Je * (Jd %*% Js)
        C <- C + a_b_one(Js, Jd) - Js * (Jd %*% Je)
        return(drop0(C + t(C)))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Id <- build_Id(adj_mat)
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        C <- a_b_one(Jd, G) - Jd * G + a_b_one(Gd, J) - Gd * J
        C <- C + Jd %*% G - Id * (Jd %*% G) + Gd %*% J - Id * (Gd %*% J)
        C <- C + a_b_one(J, Gd) - J * Gd + a_b_one(G, Jd) - G * Jd
        return(drop0(C + t(C)) / 3)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Je <- build_Je(adj_mat)
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- a_b_one(Jd, Gs) - Jd * (Gs %*% Je)
        C <- C + a_b_one(Gd, Js) - Gd * (Js %*% Je)
        C <- C + Jd %*% Gs - Je * (Jd %*% Gs) + Gd %*% Js - Je * (Gd %*% Js)
        C <- C + a_b_one(Js, Gd) - Js * (Gd %*% Je)
        C <- C + a_b_one(Gs, Jd) - Gs * (Jd %*% Je)
        return(drop0(C + t(C)) / 3)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Id <- build_Id(adj_mat)
        G <- build_G(adj_mat)
        C <- a_b_one(Gp, G) - Gp * G
        C <- C + Gp %*% G - Id * (Gp %*% G)
        C <- C + a_b_one(G, Gp) - G * Gp
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        Je <- build_Je(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- a_b_one(Gp, Gs) - Gp * (Gs %*% Je)
        C <- C + Gp %*% Gs - Je * (Gp %*% Gs)
        C <- C + a_b_one(Gs, Gp) - Gs * (Gp %*% Je)
        return(drop0(C + t(C)))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M12
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M12 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Jn <- build_Jn(adj_mat)
        J <- build_J(adj_mat)
        C <- Jd * (Jn %*% J) + Jn * (J %*% Jd) + J * (Jn %*% Jd)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        J0 <- build_J0(adj_mat)
        Js <- build_Js(adj_mat)
        C <- Jd * (J0 %*% Js) + J0 * (Js %*% Jd) + Js * (J0 %*% Jd)
        return(C + t(C))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        J <- build_J(adj_mat)
        C <- Jd * (Jn %*% G) + Gd * (Jn %*% J)
        C <- C + Jn * (J %*% Gd) + Jn * (G %*% Jd)
        C <- C + J * (Jn %*% Gd) + G * (Jn %*% Jd)
        return((C + t(C)) / 3)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        J0 <- build_J0(adj_mat)
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- Jd * (J0 %*% Gs) + Gd * (J0 %*% Js)
        C <- C + J0 * (Js %*% Gd) + J0 * (Gs %*% Jd)
        C <- C + Js * (J0 %*% Gd) + Gs * (J0 %*% Jd)
        return((C + t(C)) / 3)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- Gp * (Jn %*% G) + Jn * (G %*% Gp) + G * (Jn %*% Gp)
        return(C + t(C))
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        J0 <- build_J0(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- Gp * (J0 %*% Gs) + J0 * (Gs %*% Gp) + Gs * (J0 %*% Gp)
        return(C + t(C))
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Id <- build_Id(adj_mat)
        J <- build_J(adj_mat)
        C <- a_one_b(Jd, J) - Jd * J
        C <- C + J %*% Jd - Id * (J %*% Jd)
        C <- C + a_one_b(J, Jd) - J * Jd
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Je <- build_Je(adj_mat)
        Js <- build_Js(adj_mat)
        C <- a_one_b(Jd, Js) - Jd * (Je %*% Js)
        C <- C + Js %*% Jd - Je * (Js %*% Jd)
        C <- C + a_one_b(Js, Jd) - Js * (Je %*% Jd)
        return(drop0(C + t(C)))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Id <- build_Id(adj_mat)
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        C <- a_one_b(Jd, G) - Jd * G + a_one_b(Gd, J) - Gd * J
        C <- C + J %*% Gd - Id * (J %*% Gd) + G %*% Jd - Id * (G %*% Jd)
        C <- C + a_one_b(J, Gd) - J * Gd + a_one_b(G, Jd) - G * Jd
        return(drop0(C + t(C)) / 3)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Je <- build_Je(adj_mat)
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- a_one_b(Jd, Gs) - Jd * (Je %*% Gs)
        C <- C + a_one_b(Gd, Js) - Gd * (Je %*% Js)
        C <- C + Js %*% Gd - Je * (Js %*% Gd) + Gs %*% Jd - Je * (Gs %*% Jd)
        C <- C + a_one_b(Js, Gd) - Js * (Je %*% Gd)
        C <- C + a_one_b(Gs, Jd) - Gs * (Je %*% Jd)
        return(drop0(C + t(C)) / 3)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Id <- build_Id(adj_mat)
        G <- build_G(adj_mat)
        C <- a_one_b(Gp, G) - Gp * G
        C <- C + G %*% Gp - Id * (G %*% Gp)
        C <- C + a_one_b(G, Gp) - G * Gp
        return(drop0(C + t(C)))
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        Je <- build_Je(adj_mat)
        Gs <- build_Gs(adj_mat)
        C <- a_one_b(Gp, Gs) - Gp * (Je %*% Gs)
        C <- C + Gs %*% Gp - Je * (Gs %*% Gp)
        C <- C + a_one_b(Gs, Gp) - Gs * (Je %*% Gp)
        return(drop0(C + t(C)))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif M13
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_M13 <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jd * (Jd %*% Jn)
        Cprime <- Jn * (Jd %*% Jd)
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Jd * (Jd %*% J0)
        Cprime <- J0 * (Jd %*% Jd)
        return(C + t(C) + Cprime)
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Jn <- build_Jn(adj_mat)
        Gd <- build_Gd(adj_mat)
        C <- Jd * (Gd %*% Jn) + Gd * (Jd %*% Jn) + Jn * (Jd %*% Gd)
        return((C + t(C)) / 4)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        J0 <- build_J0(adj_mat)
        Gd <- build_Gd(adj_mat)
        C <- Jd * (Gd %*% J0) + Gd * (Jd %*% J0) + J0 * (Jd %*% Gd)
        return((C + t(C)) / 4)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Gp * (Gp %*% Jn)
        Cprime <- Jn * (Gp %*% Gp)
        return(C + t(C) + Cprime)
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- Gp * (Gp %*% J0)
        Cprime <- J0 * (Gp %*% Gp)
        return(C + t(C) + Cprime)
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(Jd, Jd) - Jd * Jd
        Cprime <- Jd %*% Jd - Id * (Jd %*% Jd)
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Jd, Jd) - Jd * (Jd %*% Je)
        Cprime <- Jd %*% Jd - Je * (Jd %*% Jd)
        return(drop0(C + t(C) + Cprime))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(Jd, Gd) - Jd * Gd + a_b_one(Gd, Jd) - Gd * Jd
        C <- C + Jd %*% Gd - Id * (Jd %*% Gd)
        return(drop0(C + t(C)) / 4)
      }

      if (motif_type == "struc") {
        Jd <- build_Jd(adj_mat)
        Gd <- build_Gd(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Jd, Gd) - Jd * (Gd %*% Je)
        C <- C + a_b_one(Gd, Jd) - Gd * (Jd %*% Je)
        C <- C + Jd %*% Gd - Je * (Jd %*% Gd)
        return(drop0(C + t(C)) / 4)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        Gp <- build_Gp(adj_mat)
        Id <- build_Id(adj_mat)
        C <- a_b_one(Gp, Gp) - Gp * Gp
        Cprime <- Gp %*% Gp - Id * (Gp %*% Gp)
        return(drop0(C + t(C) + Cprime))
      }

      if (motif_type == "struc") {
        Gp <- build_Gp(adj_mat)
        Je <- build_Je(adj_mat)
        C <- a_b_one(Gp, Gp) - Gp * (Gp %*% Je)
        Cprime <- Gp %*% Gp - Je * (Gp %*% Gp)
        return(drop0(C + t(C) + Cprime))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif Mcoll
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_Mcoll <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jn * (J %*% t(J))
        return(C)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (Js %*% t(Js))
        return(C)
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jn * (J %*% t(G)) + Jn * (G %*% t(J))
        return(C / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (Js %*% t(Gs)) + J0 * (Gs %*% t(Js))
        return(C / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jn * (G %*% t(G))
        return(C)
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (Gs %*% t(Gs))
        return(C)
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Id <- build_Id(adj_mat)
        C <- J %*% t(J) - Id * (J %*% t(J))
        return(drop0(C))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Je <- build_Je(adj_mat)
        C <- Js %*% t(Js) - Je * (Js %*% t(Js))
        return(drop0(C))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- J %*% t(G) - Id * (J %*% t(G))
        C <- C + G %*% t(J) - Id * (G %*% t(J))
        return(drop0(C) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- Js %*% t(Gs) - Je * (Js %*% t(Gs))
        C <- C + Gs %*% t(Js) - Je * (Gs %*% t(Js))
        return(drop0(C) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- G %*% t(G) - Id * (G %*% t(G))
        return(drop0(C))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- Gs %*% t(Gs) - Je * (Gs %*% t(Gs))
        return(drop0(C))
      }
    }
  }
}

#' Perform the motif adjacency matrix calculations for motif Mexpa
#'
#' @param adj_mat Adjacency matrix from which to build the motif
#' adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param mam_weight_type The weighting scheme to use.
#' One of \code{"unweighted"}, \code{"mean"} or \code{"product"}.
#' @param mam_method Which formulation to use.
#' One of \code{"dense"} or \code{"sparse"}.
#' @return A motif adjacency matrix.
#' @keywords internal
#' @importFrom Matrix drop0 t

mam_Mexpa <- function(adj_mat, motif_type, mam_weight_type, mam_method) {

  if (mam_method == "dense") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jn * (t(J) %*% J)
        return(C)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (t(Js) %*% Js)
        return(C)
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Jn <- build_Jn(adj_mat)
        G <- build_G(adj_mat)
        C <- Jn * (t(J) %*% G) + Jn * (t(G) %*% J)
        return(C / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (t(Js) %*% Gs) + J0 * (t(Gs) %*% Js)
        return(C / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Jn <- build_Jn(adj_mat)
        C <- Jn * (t(G) %*% G)
        return(C)
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        J0 <- build_J0(adj_mat)
        C <- J0 * (t(Gs) %*% Gs)
        return(C)
      }
    }
  }

  if (mam_method == "sparse") {
    if (mam_weight_type == "unweighted") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        Id <- build_Id(adj_mat)
        C <- t(J) %*% J - Id * (t(J) %*% J)
        return(drop0(C))
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Je <- build_Je(adj_mat)
        C <- t(Js) %*% Js - Je * (t(Js) %*% Js)
        return(drop0(C))
      }
    }

    if (mam_weight_type == "mean") {
      if (motif_type == "func") {
        J <- build_J(adj_mat)
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- t(J) %*% G - Id * (t(J) %*% G)
        C <- C + t(G) %*% J - Id * (t(G) %*% J)
        return(drop0(C) / 2)
      }

      if (motif_type == "struc") {
        Js <- build_Js(adj_mat)
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- t(Js) %*% Gs - Je * (t(Js) %*% Gs)
        C <- C + t(Gs) %*% Js - Je * (t(Gs) %*% Js)
        return(drop0(C) / 2)
      }
    }

    if (mam_weight_type == "product") {
      if (motif_type == "func") {
        G <- build_G(adj_mat)
        Id <- build_Id(adj_mat)
        C <- t(G) %*% G - Id * (t(G) %*% G)
        return(drop0(C))
      }

      if (motif_type == "struc") {
        Gs <- build_Gs(adj_mat)
        Je <- build_Je(adj_mat)
        C <- t(Gs) %*% Gs - Je * (t(Gs) %*% Gs)
        return(drop0(C))
      }
    }
  }
}
