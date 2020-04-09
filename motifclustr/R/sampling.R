#' Sample a DSBM
#'
#' Sample the (weighted) adjacency matrix of a (weighted) directed stochastic
#' block model (DSBM) with specified parameters.
#' @param block_sizes A vector containing the size of each block of vertices.
#' @param connection_matrix A matrix containing the block-to-block connection
#' probabilities.
#' @param weight_type The type of weighting scheme.
#' One of "unweighted", "deterministic" or "poisson".
#' @param weight_matrix A matrix containing the block-to-block weight
#' parameters.
#' Unused for weight_type="deterministic".
#' Defaults to NULL.
#' @return A randomly sampled (weighted) adjacency matrix of a DSBM.
#' @export
#' @importFrom stats rpois rbinom
#' @examples
#' block_sizes = c(10,10)
#' connection_matrix = matrix(c(0.8,0.1,0.1,0.8), nrow=2, byrow=TRUE)
#' weight_type = "poisson"
#' weight_matrix = matrix(c(10,3,3,10), nrow=2, byrow=TRUE)
#' sample_dsbm(block_sizes, connection_matrix, weight_type, weight_matrix)

sample_dsbm <- function(block_sizes, connection_matrix,
                       weight_type=c("unweighted", "deterministic", "poisson"),
                       weight_matrix=NULL){

  # check args
  if(!(is.vector(block_sizes))){
    stop("block_sizes must be a vector.")
  }
  if(!all.equal(block_sizes, as.integer(block_sizes))){
    stop("block_sizes must be integers.")
  }
  if(!all(block_sizes > 0)){
    stop("block_sizes must be at least 1.")
  }
  if(!(length(block_sizes) == nrow(connection_matrix))){
    stop("connection_matrix must have length(block_sizes) rows.")
  }
  if(!(length(block_sizes) == ncol(connection_matrix))){
    stop("connection_matrix must have length(block_sizes) columns.")
  }
  if(!(all(connection_matrix >= 0) & all(connection_matrix <= 1))){
    stop("connection_matrix entries must be in [0,1].")
  }
  weight_type <- match.arg(weight_type)
  if((weight_type != "unweighted") & is.null(weight_matrix)){
    stop("weighted requires a weight_matrix")
  }
  if(!is.null(weight_matrix)){
    if(!(length(block_sizes) == nrow(weight_matrix))){
      stop("weight_matrix must have length(block_sizes) rows.")
    }
    if(!(length(block_sizes) == ncol(weight_matrix))){
      stop("weight_matrix must have length(block_sizes) columns.")
    }
    if(!all(weight_matrix >= 0)){
      stop("weight_matrix entries must be non-negative.")
    }
  }

  # initialize variables
  n <- sum(block_sizes)
  k <- length(block_sizes)
  cumul_sizes <- c(0, cumsum(block_sizes))
  adj_mat <- matrix(0, nrow=n, ncol=n)

  for (i in 1:k) {
    for (j in 1:k) {

      # block parameters
      x_range <- (cumul_sizes[i]+1):cumul_sizes[i+1]
      y_range <- (cumul_sizes[j]+1):cumul_sizes[j+1]
      n_cells <- block_sizes[i]*block_sizes[j]
      p <- connection_matrix[i,j]

      # fill the adjacency matrix block-by-block
      adj_mat[x_range, y_range] <- rbinom(n_cells, 1, p)

      # deterministic weights
      if(weight_type == "deterministic"){
        w <- weight_matrix[i,j]
        adj_mat[x_range, y_range] <- w * adj_mat[x_range, y_range]
      }

      # poisson weights
      else if(weight_type == "poisson"){
        w <- weight_matrix[i,j]
        weights <- rpois(n_cells, w)
        adj_mat[x_range, y_range] <- weights * adj_mat[x_range, y_range]
      }
    }
  }

  # remove self-loops and make sparse
  adj_mat <- drop0_killdiag(adj_mat)

  return(adj_mat)
}

#' Sample a BSBM
#'
#' Sample the (weighted) adjacency matrix of a (weighted) bipartite stochastic
#' block model (BSBM) with specified parameters.
#' @param source_block_sizes A vector containing the size of each block
#' of source vertices.
#' @param dest_block_sizes A vector containing the size of each block
#' of destination vertices.
#' @param bipartite_connection_matrix A matrix containing the block-to-block
#' connection probabilities.
#' @param weight_type The type of weighting scheme.
#' One of "unweighted", "deterministic" or "poisson".
#' @param bipartite_weight_matrix A matrix containing the
#' sourece block-to-destination block weight parameters.
#' Unused for weight_type="deterministic".
#' Defaults to NULL.
#' @return A randomly sampled (weighted) adjacency matrix of a BSBM.
#' @export
#' @examples
#' source_block_sizes = c(10,10)
#' dest_block_sizes = c(10,10,10)
#' bipartite_connection_matrix = matrix(c(0.8,0.5,0.1,0.1,0.5,0.8),
#' nrow=2, byrow=TRUE)
#' weight_type = "poisson"
#' bipartite_weight_matrix = matrix(c(20,10,2,2,10,20), nrow=2, byrow=TRUE)
#' sample_bsbm(source_block_sizes, dest_block_sizes,
#'       bipartite_connection_matrix, weight_type, bipartite_weight_matrix)

sample_bsbm <- function(source_block_sizes, dest_block_sizes,
                       bipartite_connection_matrix,
                       weight_type=c("unweighted", "deterministic", "poisson"),
                       bipartite_weight_matrix=NULL){

  # check args
  if(!(length(source_block_sizes) == nrow(bipartite_connection_matrix))){
    stop("length(source_block_sizes) must equal
         nrow(bipartite_connection_matrix)")
  }
  if(!(length(dest_block_sizes) == ncol(bipartite_connection_matrix))){
    stop("length(dest_block_sizes) must equal
         ncol(bipartite_connection_matrix)")
  }
  if((weight_type != "unweighted") & is.null(bipartite_weight_matrix)){
    stop("weighted requires a bipartite_weight_matrix")
  }
  if(!is.null(bipartite_weight_matrix)){
    if(!(length(source_block_sizes) == nrow(bipartite_weight_matrix))){
      stop("length(source_block_sizes) must equal
           nrow(bipartite_weight_matrix)")
    }
    if(!(length(dest_block_sizes) == ncol(bipartite_weight_matrix))){
      stop("length(dest_block_sizes) must equal ncol(bipartite_weight_matrix)")
    }
  }

  # initialize parameters
  ks <- length(source_block_sizes)
  kd <- length(dest_block_sizes)
  zeros_ss <- matrix(0, nrow=ks, ncol=ks)
  zeros_d <- matrix(0, nrow=kd, ncol=(ks+kd))

  # build block sizes vector
  block_sizes <- c(source_block_sizes, dest_block_sizes)

  # build connection matrix
  connection_matrix <- cbind(zeros_ss, bipartite_connection_matrix)
  connection_matrix <- rbind(connection_matrix, zeros_d)

  # build weight matrix
  if(!is.null(bipartite_weight_matrix)){
    weight_matrix <- cbind(zeros_ss, bipartite_weight_matrix)
    weight_matrix <- rbind(weight_matrix, zeros_d)
  }
  else{
    weight_matrix <- NULL
  }

  # sample BSBM
  adj_mat <- sample_dsbm(block_sizes, connection_matrix,
                         weight_type, weight_matrix)

  return(adj_mat)
}

#' Generate a small graph for demonstrations
#'
#' Generate the sparse and dense adjacency matrices of a small weighted
#' directed graph, for demonstrating methods and running tests.
#' @return A list with two entries:
#' adj_mat_dense is the adjacency matrix in dense form.
#' adj_mat_sparse is the adjacency matrix in sparse form.

demonstration_graph<- function(){

  adj_mat_dense <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    2, 0, 3, 0, 6, 8, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 5, 0, 0, 0,14, 0, 0,18,19, 0, 0,
    0, 7, 9, 0,13, 0, 0, 0, 0, 0,21, 0,
    0, 0,11,12, 0,15, 0,17, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,16, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,24, 0, 0,
    0, 0, 0, 0, 0,20, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,22, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,23, 0, 0, 0, 0, 0
  ), nrow=12, byrow=TRUE)

  adj_mat_sparse <- drop0(adj_mat_dense)

  ans <- list(adj_mat_dense = adj_mat_dense, adj_mat_sparse = adj_mat_sparse)

  return(ans)
}
