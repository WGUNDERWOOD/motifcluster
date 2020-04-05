#' Build a motif adjacency matrix
#'
#' Build a motif adjacency matrix from an adjacency matrix.
#' @param adj_mat Adjacency matrix from which to build the motif adjacency matrix.
#' @param motif_name Motif used for the motif adjacency matrix.
#' @param motif_type Type of motif adjacency matrix to build.
#' @param weight_type The weighting scheme to use. One of "unweighted", "mean" or "product".
#' @param method Which formulation to use. One of "dense" or "sparse".
#' @return A motif adjacency matrix.
#' @keywords motif adjacency matrix
#' @importFrom Matrix drop0 t
#' @export
#' @examples

build_motif_adjacency_matrix <- function(adj_mat, motif_name, motif_type=c("func","struc"),
                                         weight_type=c("unweighted", "mean", "product"),
                                         method=c("dense", "sparse")){

  # TODO add weight type and dense/sparse method

  # check args
  if(!(motif_name %in% get_motif_names())){
    stop("Invalid motif name.")
  }
  motif_type <- match.arg(motif_type)

  if(motif_name == "Ms"){
    return(mam_Ms(adj_mat, motif_type, weight_type))
  }

  if(motif_name == "Md"){
    return(mam_Md(adj_mat, weight_type))
  }

  if(motif_name == "M1"){
    return(mam_M1(adj_mat, motif_type, weight_type))
  }

  if(motif_name == "M2"){
    return(mam_M2(adj_mat, motif_type, weight_type))
  }
}

#' Set diagonal entries to zero and sparsify
#'
#' Set the diagonal entries of a matrix to zero
#' and convert it to sparse matrix form.
#' @param mat A matrix.
#' @return A sparse-form copy of mat with its diagonal entries set to zero.
#' @importFrom Matrix drop0
#' @keywords matrix diagonal sparse

drop0_killdiag <- function(mat){

  ans <- mat
  diag(ans) <- 0
  ans <- unname(drop0(ans))

  return(ans)
}

#' Get largest connected component
#'
#' Get the indices of the largest connected component of a graph
#' from its adjacency matrix.
#' @param adj_mat An adjacency matrix of a graph.
#' @return A vector of indices corresponding to the vertices in the largest
#' connected component.

get_largest_component <- function(adj_mat){

  n <- nrow(adj_mat)
  gr <- igraph::graph_from_adjacency_matrix(adj_mat + t(adj_mat))
  comps <- igraph::components(gr)
  verts_to_keep <- (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)
}

#' Get common motif names
#'
#' Get the names of some common motifs
#' @return A vector of names of common motifs.

get_motif_names <- function(){

  motif_names = c("Ms", "Md")

  for(i in 1:13){
    motif_name = paste("M", i, sep="")
    motif_names = c(motif_names, motif_name)
  }

  motif_names = c(motif_names, c("coll", "expa", "path"))

  return(motif_names)
}
