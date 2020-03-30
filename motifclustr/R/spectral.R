computeTopSpectrum = function(L, typeLap=c('comb','rw'), topk){

  # Computes eigenvectors/values of lowest topk eigenvalues of a Laplacian L.

  if(typeLap == 'comb'){
    ansEigs = eigs_sym(L, topk, which = 'SM')
  }
  else if(typeLap == 'rw'){
    ansEigs = eigs(L, topk, which = 'SM')
  }

  inds = seq(topk,1,-1)
  vects = Re(ansEigs[['vectors']])[,inds]
  vals = Re(ansEigs[['values']])[inds]

  ansSpect = list()
  ansSpect[['vects']] = vects
  ansSpect[['vals']]  = vals

  return(ansSpect)

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
