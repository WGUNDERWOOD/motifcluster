sampleDSBM = function(block_sizes, sparsity_matrix){

  # Samples the adjacency matrix of an unweighted directed stochastic
  # block model (DSBM) as defined in the paper.

  # initialize variables
  n = sum(block_sizes)
  k = length(block_sizes)
  cumul_sizes = c(0,cumsum(block_sizes))
  G = matrix(0, nrow=n, ncol=n)

  for(i in 1:k){
    for(j in 1:k){

      # fill the matrix block-by-block
      x_range = (cumul_sizes[i]+1):cumul_sizes[i+1]
      y_range = (cumul_sizes[j]+1):cumul_sizes[j+1]
      G[x_range, y_range] = rbinom(block_sizes[i]*block_sizes[j], 1, sparsity_matrix[i,j])
    }
  }

  diag(G) = 0

  return(drop0(G))
}

sampleBSBM = function(dest_block_sizes, targ_block_sizes, p1, p2){

  # Samples the adjacency matrix of an unweighted bipartite stochastic
  # block model (BSBM) as defined in the paper.

  sparsity_matrix = matrix(c(0,0,p1,p2,
                             0,0,p2,p1,
                             0,0,0,0,
                             0,0,0,0), nrow=4, byrow=TRUE)

  G = sampleDSBM(c(dest_block_sizes, targ_block_sizes), sparsity_matrix)

  return(G)
}
