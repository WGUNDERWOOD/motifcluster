#' Run Laplace embedding
#'
#' Sample the (weighted) adjacency matrix of a (weighted) directed stochastic block
#' model (DSBM) with specified parameters.
#' @param block_sizes A vector containing the size of each block of vertices.
#' @param connection_matrix A matrix containing the block-to-block connection probabilities.
#' @param weight_type The type of weighting scheme.
#' One of \code{"unweighted"}, \code{"deterministic"} or \code{"poisson"}.
#' @param weight_matrix A matrix containing the block-to-block weight parameters.
#' Unused for \code{{weight_type="deterministic"}}.
#' Defaults to \code{NULL}.
#' @return A randomly sampled (weighted) adjacency matrix.
#' @keywords sample DSBM model
#' @export
#' @examples
# TODO example here

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
  if(!(is.matrix(connection_matrix))){
    stop("connection_matrix must be a matrix")
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
  weight_type = match.arg(weight_type)
  if((weight_type != "unweighted") & is.null(weight_matrix)){
    stop("weighted requires a weight_matrix")
  }
  if(!is.null(weight_matrix)){
    if(!(is.matrix(weight_matrix))){
      stop("weight_matrix must be a matrix")
    }
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

  for(i in 1:k){
    for(j in 1:k){

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

sampleBSBM <- function(dest_block_sizes, targ_block_sizes, p1, p2){

  # Samples the adjacency matrix of an unweighted bipartite stochastic
  # block model (BSBM) as defined in the paper.

  connection_matrix <- matrix(c(0,0,p1,p2,
                             0,0,p2,p1,
                             0,0,0,0,
                             0,0,0,0), nrow=4, byrow=TRUE)

  G <- sampleDSBM(c(dest_block_sizes, targ_block_sizes), connection_matrix)

  return(G)
}
