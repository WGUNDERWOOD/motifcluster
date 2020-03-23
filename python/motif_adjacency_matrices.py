# Functions relating to creating motif adjacency matrices

import numpy as np
import networkx as nx

def build_motif_adjacency_from_graph(G, motif_name, motif_type='struc'):
    A = nx.to_scipy_sparse_matrix(G,nodelist=sorted(G))
    return build_motif_adjacency_matrix(A, motif_name, motif_type)


def build_motif_adjacency_matrix(adjacency_matrix, motif_name, motif_type='struc'):
  # Public function
  # Builds motif adjacency matrix for a simple graph and motif.
  # motif names such as "M1" are documented in the dissertation.
  # motif_type is either: 'func'  for finding all instances of S in G.
  #                 'struc' for finding all induced instances of S in G
  if hasattr(adjacency_matrix,'todense'):
      G = adjacency_matrix.todense()
  else:
      G = adjacency_matrix

  IM = build_indicator_matrices(G)


  mac = motif_adjacency_calculations
  if(motif_type=='func'):
    motif_adjacency_matrix = mac(G, G, IM['Gd'], IM['J'], IM['Jn'], IM['J'], IM['Jd'], motif_name)
  elif(motif_type=='struc'):
      motif_adjacency_matrix = mac(G, IM['Gs'], IM['Gd'], IM['J'], IM['J0'], IM['Js'], IM['Jd'], motif_name)
#  motif_adjacency_matrix = sparse.csr_matrix(motif_adjacency_matrix)
  return motif_adjacency_matrix


def build_indicator_matrices(adjacency_matrix):
  # Builds the indicator matrices required to build a motif adjacency matrix.
  # See Section 2.2 in the dissertation for details

  G = adjacency_matrix

  J  = drop0_killdiag( (G > 0) ).astype(np.float)
  J0 = drop0_killdiag( ((G+(G.T)) == 0 ) )
  Jn = drop0_killdiag(np.array(np.ones(G.shape))) #drop0_killdiag( 1+0*G )
  Gs = drop0_killdiag( np.multiply(G,(1 - (J.T)) ))
  Gd = drop0_killdiag( np.multiply(G+G.T,np.multiply(J,(J.T))))
  Js = drop0_killdiag( (Gs > 0) ).astype(np.float)
  Jd = drop0_killdiag( (Gd > 0) ).astype(np.float)

  return {'J':J, 'J0':J0, 'Jn':Jn, 'Gs':Gs, 'Gd':Gd, 'Js':Js, 'Jd':Jd}


def motif_adjacency_calculations(G, Gs, Gd, J, J0, Js, Jd, motif_name):

  # Performs the matrix calculations which are the main part of
  # building a motif adjacency matrix.
  # See Appendix B in the dissertation for details.

  if(motif_name == 'Ms'):
    motif_adjacency_matrix = Gs + Gs.T

  elif(motif_name == 'Md'):
    motif_adjacency_matrix = Gd / 2

  elif(motif_name == 'M1'):
    C1 = np.multiply(Js.T,Js@Gs+Gs@Js)
    C2 = 0 #np.multiply(Js.T,Gs@Js)
    C3 = np.multiply(Gs.T,Js@Js)
    C = C1+C2+C3
    motif_adjacency_matrix = (C + C.T) / 3

  elif(motif_name == 'M2'):
    C  = np.multiply(Gd,Js@Js)
    C += np.multiply(Gs.T,Jd@Js + Js@Jd)
    C += np.multiply(Jd,Gs@Js + Js@Gs)
    C += np.multiply(Js.T,Gd@Js + Gs@Jd + Jd@Gs + Js@Gd)
    motif_adjacency_matrix = (C + (C.T)) / 4

  elif(motif_name == 'M3'):
    C  = np.multiply(Gd,Jd@Js + Js@Jd)
    C += np.multiply(Gs,Jd@Jd)
    C += np.multiply(Jd,Gd@Js + Gs@Jd  + Jd@Gs + Js@Gd)
    C += np.multiply(Js,Gd@Jd  + Jd@Gd)
    motif_adjacency_matrix = (C + (C.T)) / 5

  elif(motif_name == 'M4'):
    motif_adjacency_matrix = (np.multiply(Jd,Jd@Gd) + np.multiply(Jd,Gd@Jd) + np.multiply(Gd,Jd@Jd)) / 6

  elif(motif_name == 'M5'):
    C  = np.multiply(Gs,(Js.T)@Js + Js@(Js + Js.T))
    C += np.multiply(Js,(Gs.T)@Js + (Js.T)@Gs + Gs@(Js + Js.T) + Js@(Gs + Gs.T))
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M6'):
    C = np.multiply(Js,Js@Gd) + np.multiply(Js,Gs@Jd) + np.multiply(Gs,Js@Jd)
    Cprime = np.multiply(Jd,(Js.T)@Gs) + np.multiply(Jd,(Gs.T)@Js) + np.multiply(Gd,(Js.T)@Js)
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 4

  elif(motif_name == 'M7'):
    C = np.multiply(Js,Jd@Gs) + np.multiply(Js,Gd@Js) + np.multiply(Gs,Jd@Js)
    Cprime = np.multiply(Jd,Js@(Gs.T)) + np.multiply(Jd,Gs@(Js.T)) + np.multiply(Gd,Js@(Js.T))
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 4

  elif(motif_name == 'M8'):
    C = np.multiply(Js,Gs@J0) + np.multiply(Gs,Js@J0)
    Cprime = np.multiply(J0,(Js.T)@Gs) + np.multiply(J0,(Gs.T)@Js)
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 2

  elif(motif_name == 'M9'):
    C = np.multiply(Js,J0@Gs.T) + np.multiply(Gs,J0@(Js.T))
    C = C + np.multiply(J0,Js@Gs) + np.multiply(J0,Gs@Js)
    C = C + np.multiply(Js,(Gs.T)@J0) + np.multiply(Gs,(Js.T)@J0)
    motif_adjacency_matrix = (C + (C.T)) / 2

  elif(motif_name == 'M10'):
    C = np.multiply(Js,J0@Gs) + np.multiply(Gs,J0@Js)
    Cprime = np.multiply(J0,Js@(Gs.T)) + np.multiply(J0,Gs@(Js.T))
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 2

  elif(motif_name == 'M11'):
    C = np.multiply(Jd,Gs@J0) + np.multiply(Gd,Js@J0)
    C = C + np.multiply(J0,Jd@Gs) + np.multiply(J0,Gd@Js)
    C = C + np.multiply(Js,Gd@J0) + np.multiply(Gs,Jd@J0)
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M12'):
    C = np.multiply(Jd,J0@Gs) + np.multiply(Gd,J0@Js)
    C = C + np.multiply(J0,Js@Gd) + np.multiply(J0,Gs@Jd)
    C = C + np.multiply(Js,J0@Gd) + np.multiply(Gs,J0@Jd)
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M13'):
    C = np.multiply(Jd,Gd@J0) + np.multiply(Gd,Jd@J0) + np.multiply(J0,Jd@Gd)
    motif_adjacency_matrix = (C + (C.T)) / 4

  elif(motif_name == 'coll'):
    C = np.multiply(J0,Js@(Gs.T)) + np.multiply(J0,Gs@(Js.T))
    motif_adjacency_matrix = C / 2

  elif(motif_name == 'expa'):
    C = np.multiply(J0,(Js.T)@Gs) + np.multiply(J0,(Gs.T)@Js)
    motif_adjacency_matrix = C / 2

  elif(motif_name == 'path'):
    C = np.multiply(J0,Js@Gs) + np.multiply(J0,Gs@Js)
    motif_adjacency_matrix = (C + (C.T)) / 2


  return(motif_adjacency_matrix)


def drop0_killdiag(some_matrix):
    M = some_matrix
    for i in range(M.shape[0]):
        M[i,i] = 0
    return M
