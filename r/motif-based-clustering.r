# Generative Models
sampleDSBM = function(block_sizes, sparsity_matrix){

  n = sum(block_sizes)
  k = length(block_sizes)
  cumul_sizes = c(0,cumsum(block_sizes))
  G = matrix(0, nrow=n, ncol=n)

  for(i in 1:k){
    for(j in 1:k){
      x_range = (cumul_sizes[i]+1):cumul_sizes[i+1]
      y_range = (cumul_sizes[j]+1):cumul_sizes[j+1]
      G[x_range, y_range] = rbinom(block_sizes[i]*block_sizes[j], 1, sparsity_matrix[i,j])
    }
  }

  diag(G) = 0

  return(drop0(G))
}

sampleBSBM = function(dest_block_sizes, targ_block_sizes, p1, p2){

  sparsity_matrix = matrix(c(0,0,p1,p2,
                             0,0,p2,p1,
                             0,0,0,0,
                             0,0,0,0), nrow=4, byrow=TRUE)

  G = sampleDSBM(c(dest_block_sizes, targ_block_sizes), sparsity_matrix)

  return(G)
}


# Motif adjacency matrices
motifAdjacency = function(G, motif_name, type=c('func','struc')){

  # Builds motif adjacency matrix for a simple graph adjacency matrix G and motif name.
  # type is either: 'func'  for finding all instances of S in G.
  #                 'struc' for finding all induced instances of S in G

  type = match.arg(type)
  IM = buildIndMats(G)

  if(type=='func'){

    motifadj = motifAdjCalcs(G, G, IM$Gd, IM$J, IM$Jn, IM$J, IM$Jd, motif_name)
  }

  else if(type=='struc'){

    motifadj = motifAdjCalcs(G, IM$Gs, IM$Gd, IM$J, IM$J0, IM$Js, IM$Jd, motif_name)
  }

  motifadj = unname(drop0(motifadj))

  return(motifadj)
}

buildIndMats = function(G){

  J  = drop0_killdiag( 1*(G > 0) )
  J0 = drop0_killdiag( 1*((G+t(G)) == 0 ) )
  Jn = drop0_killdiag( 1+0*G )
  Gs = drop0_killdiag( G*(1 - t(J)) )
  Gd = drop0_killdiag( (G+t(G)) * J * t(J) )
  Js = drop0_killdiag( 1*(Gs > 0) )
  Jd = drop0_killdiag( 1*(Gd > 0) )

  return(list(J=J, J0=J0, Jn=Jn, Gs=Gs, Gd=Gd, Js=Js, Jd=Jd))
}

motifAdjCalcs = function(G, Gs, Gd, J, J0, Js, Jd, motif_name){

  if(motif_name == 'Ms'){
    motifadj = Gs + t(Gs)
  }

  else if(motif_name == 'Md'){
    motifadj = Gd / 2
  }

  else if(motif_name == 'M1'){
    C = t(Js)*(Js%*%Gs) + t(Js)*(Gs%*%Js) + t(Gs)*(Js%*%Js)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M2'){
    C = t(Js)*(Jd%*%Gs) + t(Js)*(Gd%*%Js) + t(Gs)*(Jd%*%Js)
    C = C + t(Js)*(Js%*%Gd) + t(Js)*(Gs%*%Jd) + t(Gs)*(Js%*%Jd)
    C = C + Jd*(Js%*%Gs) + Jd*(Gs%*%Js) + Gd*(Js%*%Js)
    motifadj = (C + t(C)) / 4
  }

  else if(motif_name == 'M3'){
    C = Js*(Jd%*%Gd) + Js*(Gd%*%Jd) + Gs*(Jd%*%Jd)
    C = C + Jd*(Jd%*%Gs) + Jd*(Gd%*%Js) + Gd*(Jd%*%Js)
    C = C + Jd*(Js%*%Gd) + Jd*(Gs%*%Jd) + Gd*(Js%*%Jd)
    motifadj = (C + t(C)) / 5
  }

  else if(motif_name == 'M4'){
    motifadj = (Jd*(Jd%*%Gd) + Jd*(Gd%*%Jd) + Gd*(Jd%*%Jd)) / 6
  }

  else if(motif_name == 'M5'){
    C = Js*(Js%*%Gs) + Js*(Gs%*%Js) + Gs*(Js%*%Js)
    C = C + Js*(Js%*%t(Gs)) + Js*(Gs%*%t(Js)) + Gs*(Js%*%t(Js))
    C = C + Js*(t(Js)%*%Gs) + Js*(t(Gs)%*%Js) + Gs*(t(Js)%*%Js)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M6'){
    C = Js*(Js%*%Gd) + Js*(Gs%*%Jd) + Gs*(Js%*%Jd)
    Cprime = Jd*(t(Js)%*%Gs) + Jd*(t(Gs)%*%Js) + Gd*(t(Js)%*%Js)
    motifadj = (C + t(C) + Cprime) / 4
  }

  else if(motif_name == 'M7'){
    C = Js*(Jd%*%Gs) + Js*(Gd%*%Js) + Gs*(Jd%*%Js)
    Cprime = Jd*(Js%*%t(Gs)) + Jd*(Gs%*%t(Js)) + Gd*(Js%*%t(Js))
    motifadj = (C + t(C) + Cprime) / 4
  }

  else if(motif_name == 'M8'){
    C = Js*(Gs%*%J0) + Gs*(Js%*%J0)
    Cprime = J0*(t(Js)%*%Gs) + J0*(t(Gs)%*%Js)
    motifadj = (C + t(C) + Cprime) / 2
  }

  else if(motif_name == 'M9'){
    C = Js*(J0%*%t(Gs)) + Gs*(J0%*%t(Js))
    C = C + J0*(Js%*%Gs) + J0*(Gs%*%Js)
    C = C + Js*(t(Gs)%*%J0) + Gs*(t(Js)%*%J0)
    motifadj = (C + t(C)) / 2
  }

  else if(motif_name == 'M10'){
    C = Js*(J0%*%Gs) + Gs*(J0%*%Js)
    Cprime = J0*(Js%*%t(Gs)) + J0*(Gs%*%t(Js))
    motifadj = (C + t(C) + Cprime) / 2
  }

  else if(motif_name == 'M11'){
    C = Jd*(Gs%*%J0) + Gd*(Js%*%J0)
    C = C + J0*(Jd%*%Gs) + J0*(Gd%*%Js)
    C = C + Js*(Gd%*%J0) + Gs*(Jd%*%J0)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M12'){
    C = Jd*(J0%*%Gs) + Gd*(J0%*%Js)
    C = C + J0*(Js%*%Gd) + J0*(Gs%*%Jd)
    C = C + Js*(J0%*%Gd) + Gs*(J0%*%Jd)
    motifadj = (C + t(C)) / 3
  }

  else if(motif_name == 'M13'){
    C = Jd*(Gd%*%J0) + Gd*(Jd%*%J0) + J0*(Jd%*%Gd)
    motifadj = (C + t(C)) / 4
  }

  else if(motif_name == 'coll'){
    C = J0*(Js%*%t(Gs)) + J0*(Gs%*%t(Js))
    motifadj = C / 2
  }

  else if(motif_name == 'expa'){
    C = J0*(t(Js)%*%Gs) + J0*(t(Gs)%*%Js)
    motifadj = C / 2
  }

  else if(motif_name == 'path'){
    C = J0*(Js%*%Gs) + J0*(Gs%*%Js)
    motifadj = (C + t(C)) / 2
  }


  return(motifadj)

}

drop0_killdiag = function(mat){

  ans = mat
  diag(ans) = 0
  ans = drop0(ans)

  return(ans)
}

largestComponent = function(G){

  n = nrow(G)
  Gr = graph_from_adjacency_matrix(G, weighted=TRUE)
  comps = components(Gr)
  verts_to_keep = (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)

}


# Spectral methods
computeTopSpectrum = function(L, typeLap=c('comb','rw'), topk){

  # Computes eigenvectors/values of lowest 'topk' eigenvalues of a Laplacian L.
  # Returns a list of $vects and $vals

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

  # Builds various types of Laplacian Matrix, given an adjacency matrix G and a type of Laplacian.

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
