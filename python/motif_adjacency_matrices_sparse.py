import networkx as nx
from scipy import sparse
import numpy as np


# Functions relating to creating motif adjacency matrices
def build_motif_adjacency_from_graph(G, motif_name, motif_type='struc'):
    A = nx.to_scipy_sparse_matrix(G,nodelist=sorted(G))
    return build_motif_adjacency_matrix(A, motif_name, motif_type)


def build_motif_adjacency_matrix(adjacency_matrix, motif_name, motif_type='struc'):
  # Public function
  # Builds motif adjacency matrix for a simple graph and motif.
  # motif names such as "M1" are documented in the dissertation.
  # motif_type is either: 'func'  for finding all instances of S in G.
  #                 'struc' for finding all induced instances of S in G

  G = adjacency_matrix

  IM = build_indicator_matrices(G)


  mac = motif_adjacency_calculations
  if(motif_type=='func'):
    motif_adjacency_matrix =   mac(G , G        , IM['Gd'] , IM['J'] , IM['J']  , IM['Jd'] , motif_name , motif_type)
  elif(motif_type=='struc'):
      motif_adjacency_matrix = mac(G , IM['Gs'] , IM['Gd'] , IM['J'] , IM['Js'] , IM['Jd'] , motif_name , motif_type)
#  motif_adjacency_matrix = sparse.csr_matrix(motif_adjacency_matrix)
  return motif_adjacency_matrix

def build_indicator_matrices(adjacency_matrix):
  # Builds the indicator matrices required to build a motif adjacency matrix.
  # See Section 2.2 in the dissertation for details
  G = adjacency_matrix

  I = sparse.identity(G.shape[0])
  J  = ( (G > 0) ).astype(np.float)
  J = J - I.multiply(J)
  Gs = (G - G.multiply(J.T))
  Gs = Gs - I.multiply(Gs) # This kills the diagonal
  temp = G+G.T
  temp1 = J.multiply(J.T)
  Gd = ( temp1.multiply(temp))
  Gd = Gd - I.multiply(Gd) # This kills the diagonal
  Js = (Gs > 0).astype(np.float)
  Js = Js - I.multiply(Js) # This kills the diagonal
  Jd = (Gd > 0).astype(np.float)
  Jd = Jd - I.multiply(Jd) # This kills the diagonal
  return {'J':J, 'Gs':Gs, 'Gd':Gd, 'Js':Js, 'Jd':Jd}


def motif_adjacency_calculations(G, Gs, Gd, J, Js, Jd, motif_name,motif_type=None):

  # Performs the matrix calculations which are the main part of
  # building a motif adjacency matrix.
  # See Appendix B in the dissertation for details.

  if(motif_name == 'Ms'):
    motif_adjacency_matrix = Gs + Gs.T

  elif(motif_name == 'Md'):
    motif_adjacency_matrix = Gd / 2

  elif(motif_name == 'M1'):
    C = (Js.T).multiply(Gs@Js) +  (Js.T).multiply(Js@Gs) + (Gs.T).multiply(Js@Js)
    motif_adjacency_matrix = (C + C.T) / 3

  elif(motif_name == 'M2'):
    C =     (Js.T).multiply(Jd@Gs+Js@Gd+Gs@Jd + Gd@Js)
    C = C + (Gs.T).multiply(Js@Jd + Jd@Js)
    C = C + (Jd).multiply(Js@Gs + Gs@Js) + Gd.multiply(Js@Js)
    motif_adjacency_matrix = (C + (C.T)) / 4

  elif(motif_name == 'M3'):
    C1 = Js.multiply(Jd@Gd + Gd@Jd) + Gs.multiply(Jd@Jd)
    C2 = Jd.multiply(Jd@Gs + Gd@Js + Js@Gd + Gs@Jd)
    C3 = Gd.multiply(Jd@Js + Js@Jd)
    C = C1 + C2 + C3
    motif_adjacency_matrix = (C + (C.T)) / 5

  elif(motif_name == 'M4'):
    motif_adjacency_matrix = (Jd.multiply(Jd@Gd + Gd@Jd) + Gd.multiply(Jd@Jd)) / 6

  elif(motif_name == 'M5'):
    C1 =  Js.multiply(Gs@(Js+Js.T) + Js@(Gs.T) + (Js + Js.T)@Gs + (Gs.T)@Js)
    C2 =  Gs.multiply(Js@(Js + Js.T) + (Js.T)@Js)
    C3 = 0
    C = C1 + C2 + C3
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M6'):
    C = Js.multiply(Js@Gd + Gs@Jd) + Gs.multiply(Js@Jd)
    Cprime = Jd.multiply((Js.T)@Gs + (Gs.T)@Js) + Gd.multiply((Js.T)@Js)
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 4

  elif(motif_name == 'M7'):
    C = Js.multiply(Jd@Gs + Gd@Js) + Gs.multiply(Jd@Js)
    Cprime = Jd.multiply(Js@(Gs.T) + Gs@(Js.T)) + Gd.multiply(Js@(Js.T))
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 4

  elif(motif_name == 'M8'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        C = matmultiplyBy1Right(Js,Gs,o1)+ matmultiplyBy1Right(Gs,Js,o1) -2*Gs
        t1 = (Js.T)@Gs + (Gs.T)@Js
        Cprime = t1 - I.multiply(t1)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        C = matmultiplyBy1Right(Js,Gs,o1)+ matmultiplyBy1Right(Gs,Js,o1)- Js.multiply(Gs@J1) - Gs.multiply(Js@J1)
        t1 = (Js.T)@Gs + (Gs.T)@Js
        Cprime = t1 - J1.multiply(t1)
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 2

  elif(motif_name == 'M9'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        C1 =  matmultiplyBy1Left(Js,Gs.T,o1)  + matmultiplyBy1Left(Gs,Js.T,o1)-2*Gs.T.multiply(Js)
        C2 =  matmultiplyBy1Right(Js,Gs.T,o1) + matmultiplyBy1Right(Gs,Js.T,o1) - 2*Gs.multiply(Js.T)
        t5 = Js@Gs + Gs@Js
        C3 =  t5 - I.multiply(t5)
        C = C1 + C2 + C3
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        C1  =   matmultiplyBy1Left(Js,Gs.T,o1)+ matmultiplyBy1Left(Gs,Js.T,o1) - Js.multiply(J1@Gs.T + Gs.T@J1)
        C2 =   matmultiplyBy1Right(Js,Gs.T,o1)+ matmultiplyBy1Right(Gs,Js.T,o1) - Gs.multiply(J1@Js.T + Js.T@J1)
        t5 = Js@Gs + Gs@Js
        C3 =  t5 - J1.multiply(t5)
        C = C1 + C2 + C3
    motif_adjacency_matrix = (C + (C.T)) / 2

  elif(motif_name == 'M10'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        C =      matmultiplyBy1Left(Js,Gs,o1)+ matmultiplyBy1Left(Gs,Js,o1) - Js.multiply(Gs) - Gs.multiply(Js)
        t3 =Js@(Gs.T)+ Gs@(Js.T)
        Cprime = t3 - I.multiply(t3)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        C =      matmultiplyBy1Left(Js,Gs,o1)+ matmultiplyBy1Left(Gs,Js,o1) - Js.multiply(J1@Gs) - Gs.multiply(J1@Js)
        t3 =Js@(Gs.T)+ Gs@(Js.T)
        Cprime = t3 - J1.multiply(t3)
    motif_adjacency_matrix = (C + (C.T) + Cprime) / 2

  elif(motif_name == 'M11'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        C1 =  matmultiplyBy1Right(Jd,Gs,o1)+ matmultiplyBy1Right(Gd,Js,o1) - Jd.multiply(Gs) - Gd.multiply(Js)
        t1 = Jd@Gs + Gd@Js
        C2 = t1-I.multiply(t1)
        C3 = matmultiplyBy1Right(Js,Gd,o1)+ matmultiplyBy1Right(Gs,Jd,o1) - Js.multiply(Gd) - Gs.multiply(Jd)
        C = C1 + C2 + C3
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        C1 =  matmultiplyBy1Right(Jd,Gs,o1)+ matmultiplyBy1Right(Gd,Js,o1) - Jd.multiply(Gs@J1) - Gd.multiply(Js@J1)
        # Trick to the remove elements could also do t3[J1]=0 latter might be
        # faster
        t3 = Jd@Gs + Gd@Js
        C2 = t3-J1.multiply(t3)
        C3 = matmultiplyBy1Right(Js,Gd,o1)+ matmultiplyBy1Right(Gs,Jd,o1) - Js.multiply(Gd@J1) - Gs.multiply(Jd@J1)
        C = C1 + C2 + C3
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M12'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        C =     matmultiplyBy1Left(Jd,Gs,o1)+ matmultiplyBy1Left(Gd,Js,o1) - Jd.multiply(Gs) - Gd.multiply(Js)
        t1 = Js@Gd + Gs@Jd
        C = C + t1 - I.multiply(t1)
        C = C + matmultiplyBy1Left(Js,Gd,o1) + matmultiplyBy1Left(Gs,Jd,o1) - Js.multiply(Gd) - Gs.multiply(Jd)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        C =  matmultiplyBy1Left(Jd,Gs,o1)+ matmultiplyBy1Left(Gd,Js,o1) - Jd.multiply(J1@Gs) - Gd.multiply(J1@Js)
        t1 = Js@Gd + Gs@Jd
        C = C + t1 - J1.multiply(t1)
        C = C + matmultiplyBy1Left(Js,Gd,o1) + matmultiplyBy1Left(Gs,Jd,o1) - Js.multiply(J1@Gd) - Gs.multiply(J1@Jd)
    motif_adjacency_matrix = (C + (C.T)) / 3

  elif(motif_name == 'M13'):
    I = sparse.identity(Js.shape[0])
    o1 = np.ones(I.shape[0])
    if motif_type=='func':
        t1 = Jd@Gd
        C = matmultiplyBy1Right(Jd,Gd,o1) + matmultiplyBy1Right(Gd,Jd,o1) + t1 - I.multiply(t1) - Jd.multiply(Gd) - Gd.multiply(Jd)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        t1 = Jd@Gd
        C1 = t1 - J1.multiply(t1)
        C2 = matmultiplyBy1Right(Jd,Gd,o1) + matmultiplyBy1Right(Gd,Jd,o1) - Jd.multiply(Gd@J1) - Gd.multiply(Jd@J1)
        C = C1 + C2
    motif_adjacency_matrix = (C + (C.T)) / 4

  elif(motif_name == 'coll'):
    I = sparse.identity(Js.shape[0])
    if motif_type=='func':
        t1 = Js@(Gs.T) + Gs@(Js.T)
        C = t1 - I.multiply(t1)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        t1 = Js@(Gs.T) + Gs@(Js.T)
        C = t1 - J1.multiply(t1)
    motif_adjacency_matrix = C / 2

  elif(motif_name == 'expa'):
    I = sparse.identity(Js.shape[0])
    if motif_type=='func':
        t1 = (Js.T)@Gs + (Gs.T)@Js
        C = t1 - I.multiply(t1)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        t1 = (Js.T)@Gs + (Gs.T)@Js
        C = t1 - J1.multiply(t1)
    motif_adjacency_matrix = C / 2

  elif(motif_name == 'path'):
    I = sparse.identity(Js.shape[0])
    if motif_type=='func':
        t1 = Js@Gs + Gs@Js
        C = t1 - I.multiply(t1)
    elif motif_type=='struc':
        J1 = ((I+J+J.T)>0).astype(np.float)
        t1 = Js@Gs + Gs@Js
        C = t1 - J1.multiply(t1)
    motif_adjacency_matrix = (C + (C.T)) / 2

  return(motif_adjacency_matrix)

def matmultiplyBy1Left(X,A,oneVec=None):
    if oneVec is None:
        oneVec = np.ones(X.shape[0])
    return X.multiply(oneVec@A)

def matmultiplyBy1Right(X,A,oneVec=None):
    if oneVec is None:
        oneVec = np.ones(X.shape[0])
    t1 = (X.T.multiply((A@oneVec))).T
    return t1





