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

  # Builds the indicator and variants of
  # adjacency matrices for a graph adjacency matrix G

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

  # Performs the matrix operations required to build
  # a motif adjacency matrix.

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

  # Set the diagonal entries of a matrix to zero,
  # and convert it to sparse form.

  ans = mat
  diag(ans) = 0
  ans = drop0(ans)

  return(ans)
}

largestComponent = function(G){

  # Return the vertices in the largest connected component
  # associated with an adjacency matrix.

  n = nrow(G)
  Gr = graph_from_adjacency_matrix(G, weighted=TRUE)
  comps = components(Gr)
  verts_to_keep = (1:n)[comps$membership == which.max(comps$csize)]

  return(verts_to_keep)
}
