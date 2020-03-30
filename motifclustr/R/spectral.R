#' Compute first few eigenvalues and eigenvectors.
#'
#' Compute the first few eigenvalues (by magnitude) and
#' associated eigenvectors of a matrix.
#' @param mat The symmetric matrix for which eigenvalues and eigenvectors are to be calculated.
#' @param l The number of eigenvalues and eigenvectors to calculate.
#' @return The first \code{l} eigenvalues (by magnitude) and associated eigenvectors of \code{mat}.
#' @return A list with two entries: \code{vals} contains the a vector of the first few eigenvalues,
#' and \code{vects} contains a matrix of the associated eigenvectors.
#' @keywords eigenvalue eigenvector spectrum matrix
#' @export
#' @examples
#' get_first_eigs(matrix(rep(1,9), nrow=3), 2)

get_first_eigs <- function(mat, l){

  # check args
  if(!all.equal(l, as.integer(l))){
    stop("l must be an integer.")
  }
  if(!(l > 0)){
    stop("l must be at least 1.")
  }
  if(!is.matrix(mat)){
    stop("mat must be a matrix.")}
  if(!isSymmetric(mat)){
    stop("mat must be symmetric.")
  }

  # get spectrum
  ans_eigs <- eigs(mat, l, which = 'SM')

  # order eigenvalues and eigenvectors
  inds <- seq(l, 1, -1)
  vects <- Re(ansEigs[['vectors']])[,inds]
  vals <- Re(ansEigs[['values']])[inds]

  # return a list
  ans_spect <- list()
  ans_spect[['vects']] = vects
  ans_spect[['vals']]  = vals

  return(ans_spect)
}

buildLaplacian = function(G, typeLap=c('comb','rw')){

  # Builds various types of Laplacian matrix,
  # given an adjacency matrix G and a type of Laplacian.

  # typeLap
  # 'comb'      combinatorial Laplacian                           L = D - G
  # 'rw'        random walk Laplacian (row-normalised)            L = D^(-1) (D - G)

  typeLap = match.arg(typeLap)

  degsG = apply(G, 1, sum)
  n = nrow(G)

  if (typeLap=='comb'){
    D = diag(degsG)
    L =  D - G
  }
  else if (typeLap == 'rw'){
    Dinv = diag(degsG^(-1))
    L =  diag(n) - Dinv %*% G
  }

  return(L)
}

runLapEmb = function(G, topk, typeLap=c('comb','rw')){

  # Runs Laplace embedding of an adjacency matrix G, given number of clusters and
  # Laplacian type. Returns a list with $vects and $vals.

  typeLap = match.arg(typeLap)

  L = buildLaplacian(G, typeLap)
  ansSpect = computeTopSpectrum(L, typeLap, topk)

  return(ansSpect)

}

runMotifEmb = function(G, motif_name, type, typeLap, numEigs){

  # Runs motif-based embedding on a given adjacency matrix.

  M = motifAdjacency(G, motif_name, type)
  comps = largestComponent(M)
  M = M[comps,comps, drop=FALSE]
  Gcomps = G[comps,comps, drop=FALSE]
  spect = runLapEmb(M,numEigs,typeLap)
  ans = list()
  ans$vects = spect$vects
  ans$vals = spect$vals
  ans$comps = comps
  ans$M = M
  ans$Gcomps = Gcomps

  return(ans)
}
