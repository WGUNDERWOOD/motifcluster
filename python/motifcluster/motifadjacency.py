"""
Functions for building motif adjacency matrices
are in `motifcluster.motifadjacency`.
"""

from scipy import sparse

from motifcluster import utils as mcut
from motifcluster import indicators as mcin

def build_motif_adjacency_matrix(adj_mat, motif_name, motif_type="struc",
                                 mam_weight_type="unweighted", mam_method="sparse"):

  """
  Build a motif adjacency matrix.

  Build a motif adjacency matrix from an adjacency matrix.
  Entry (`i, j`) of a motif adjacency matrix is the
  sum of the weights of all motifs containing both
  nodes `i` and `j`.

  - The motif is specified by name and the type of motif instance can be one of:

    - Functional: motifs should appear as subgraphs.
    - Structural: motifs should appear as induced subgraphs.

  - The weighting scheme can be one of:

    - Unweighted: the weight of any motif instance is one.
    - Mean: the weight of any motif instance
      is the mean of its edge weights.
    - Product: the weight of any motif instance
      is the product of its edge weights.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_name : str
    Motif used for the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
    One of `"func"` or `"struc"`.
  mam_weight_type : str
    The weighting scheme to use.
    One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use.
    One of `"dense"` or `"sparse"`.
    The sparse formulation avoids generating large dense matrices
    so tends to be faster for large sparse graphs.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.

  Examples
  --------
  >>> adj_mat = np.array(range(1, 10)).reshape((3, 3))
  >>> build_motif_adjacency_matrix(adj_mat, "M1", "func", "mean")
  """

  # check args
  assert motif_name in mcut.get_motif_names()
  assert motif_type in ["struc", "func"]
  assert mam_weight_type in ["unweighted", "mean", "product"]
  assert mam_method in ["sparse", "dense"]

  # ensure correct sparsity type
  if mam_method == "dense":
    if sparse.issparse(adj_mat):
      adj_mat = adj_mat.toarray()

  elif mam_method == "sparse":
    adj_mat = sparse.csr_matrix(adj_mat)

  # build mam
  if motif_name == "Ms":
    motif_adj_mat = mam_Ms(adj_mat, motif_type, mam_weight_type)

  if motif_name == "Md":
    motif_adj_mat = mam_Md(adj_mat, mam_weight_type)

  if motif_name == "M1":
    motif_adj_mat = mam_M1(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M2":
    motif_adj_mat = mam_M2(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M3":
    motif_adj_mat = mam_M3(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M4":
    motif_adj_mat = mam_M4(adj_mat, mam_weight_type, mam_method)

  if motif_name == "M5":
    motif_adj_mat = mam_M5(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M6":
    motif_adj_mat = mam_M6(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M7":
    motif_adj_mat = mam_M7(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M8":
    motif_adj_mat = mam_M8(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M9":
    motif_adj_mat = mam_M9(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M10":
    motif_adj_mat = mam_M10(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M11":
    motif_adj_mat = mam_M11(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M12":
    motif_adj_mat = mam_M12(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "M13":
    motif_adj_mat = mam_M13(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "Mcoll":
    motif_adj_mat = mam_Mcoll(adj_mat, motif_type, mam_weight_type, mam_method)

  if motif_name == "Mexpa":
    motif_adj_mat = mam_Mexpa(adj_mat, motif_type, mam_weight_type, mam_method)

  return sparse.csr_matrix(motif_adj_mat)


def mam_Ms(adj_mat, motif_type, mam_weight_type):

  """
  Perform the motif adjacency matrix calculations for motif Ms.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_weight_type == "unweighted":
    if motif_type == "func":
      J = mcin._build_J(adj_mat)
      motif_adj_mat = J + J.transpose()

    if motif_type == "struc":
      Js = mcin._build_Js(adj_mat)
      motif_adj_mat = Js + Js.transpose()

  if mam_weight_type == "mean":
    if motif_type == "func":
      G = mcin._build_G(adj_mat)
      motif_adj_mat = G + G.transpose()

    if motif_type == "struc":
      Gs = mcin._build_Gs(adj_mat)
      motif_adj_mat = Gs + Gs.transpose()

  if mam_weight_type == "product":
    if motif_type == "func":
      G = mcin._build_G(adj_mat)
      motif_adj_mat = G + G.transpose()

    if motif_type == "struc":
      Gs = mcin._build_Gs(adj_mat)
      motif_adj_mat = Gs + Gs.transpose()

  return motif_adj_mat


def mam_Md(adj_mat, mam_weight_type):

  """
  Perform the motif adjacency matrix calculations for motif Md.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_weight_type == "unweighted":
    Jd = mcin._build_Jd(adj_mat)
    motif_adj_mat = Jd

  if mam_weight_type == "mean":
    Gd = mcin._build_Gd(adj_mat)
    motif_adj_mat = Gd / 2

  if mam_weight_type == "product":
    Gp = mcin._build_Gp(adj_mat)
    motif_adj_mat = Gp

  return motif_adj_mat


def mam_M1(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M1.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        C = J.transpose() * (J @ J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        C = Js.transpose() * (Js @ Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.transpose() * (J @ G) + J.transpose() * (G @ J)
        C += G.transpose() * (J @ J)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Js.transpose() * (Js @ Gs) + Js.transpose() * (Gs @ Js)
        C += Gs.transpose() * (Js @ Js)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        C = G.transpose() * (G @ G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        C = Gs.transpose() * (Gs @ Gs)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        C = J.transpose().multiply(J * J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        C = Js.transpose().multiply(Js * Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.transpose().multiply(J * G) + J.transpose().multiply(G * J)
        C += G.transpose().multiply(J * J)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Js.transpose().multiply(Js * Gs) + Js.transpose().multiply(Gs * Js)
        C += Gs.transpose().multiply(Js * Js)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        C = G.transpose().multiply(G * G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        C = Gs.transpose().multiply(Gs * Gs)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M2(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements

  """
  Perform the motif adjacency matrix calculations for motif M2.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J.transpose() * (Jd @ J) + J.transpose() * (J @ Jd)
        C += Jd * (J @ J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js.transpose() * (Jd @ Js) + Js.transpose() * (Js @ Jd)
        C += Jd * (Js @ Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.transpose() * (Jd @ G) + J.transpose() * (Gd @ J)
        C += G.transpose() * (Jd @ J)
        C = C + J.transpose() * (J @ Gd) + J.transpose() * (G @ Jd)
        C += G.transpose() * (J @ Jd)
        C = C + Jd * (J @ G) + Jd * (G @ J) + Gd * (J @ J)
        motif_adj_mat = (C + C.transpose()) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js.transpose() * (Jd @ Gs) + Js.transpose() * (Gd @ Js)
        C += Gs.transpose() * (Jd @ Js)
        C = C + Js.transpose() * (Js @ Gd)
        C += Js.transpose() * (Gs @ Jd) + Gs.transpose() * (Js @ Jd)
        C = C + Jd * (Js @ Gs) + Jd * (Gs @ Js) + Gd * (Js @ Js)
        motif_adj_mat = (C + C.transpose()) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G.transpose() * (Gp @ G) + G.transpose() * (G @ Gp)
        C += Gp * (G @ G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs.transpose() * (Gp @ Gs) + Gs.transpose() * (Gs @ Gp)
        C += Gp * (Gs @ Gs)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J.transpose().multiply(Jd * J) + J.transpose().multiply(J * Jd)
        C += Jd.multiply(J * J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js.transpose().multiply(Jd * Js) + Js.transpose().multiply(Js * Jd)
        C += Jd.multiply(Js * Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.transpose().multiply(Jd * G) + J.transpose().multiply(Gd * J)
        C += G.transpose().multiply(Jd * J)
        C = C + J.transpose().multiply(J * Gd) + J.transpose().multiply(G * Jd)
        C += G.transpose().multiply(J * Jd)
        C = C + Jd.multiply(J * G) + Jd.multiply(G * J) + Gd.multiply(J * J)
        motif_adj_mat = (C + C.transpose()) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js.transpose().multiply(Jd * Gs) + Js.transpose().multiply(Gd * Js)
        C += Gs.transpose().multiply(Jd * Js)
        C = C + Js.transpose().multiply(Js * Gd)
        C += Js.transpose().multiply(Gs * Jd) + Gs.transpose().multiply(Js * Jd)
        C = C + Jd.multiply(Js * Gs) + Jd.multiply(Gs * Js) + Gd.multiply(Js * Js)
        motif_adj_mat = (C + C.transpose()) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G.transpose().multiply(Gp * G) + G.transpose().multiply(G * Gp)
        C += Gp.multiply(G * G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs.transpose().multiply(Gp * Gs) + Gs.transpose().multiply(Gs * Gp)
        C += Gp.multiply(Gs * Gs)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M3(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M3.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J * (Jd @ Jd) + Jd * (Jd @ J) + Jd * (J @ Jd)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js * (Jd @ Jd) + Jd * (Jd @ Js) + Jd * (Js @ Jd)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (Jd @ Gd) + J * (Gd @ Jd) + G * (Jd @ Jd)
        C = C + Jd * (Jd @ G) + Jd * (Gd @ J) + Gd * (Jd @ J)
        C = C + Jd * (J @ Gd) + Jd * (G @ Jd) + Gd * (J @ Jd)
        motif_adj_mat = (C + C.transpose()) / 5

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js * (Jd @ Gd) + Js * (Gd @ Jd) + Gs * (Jd @ Jd)
        C = C + Jd * (Jd @ Gs) + Jd * (Gd @ Js) + Gd * (Jd @ Js)
        C = C + Jd * (Js @ Gd) + Jd * (Gs @ Jd) + Gd * (Js @ Jd)
        motif_adj_mat = (C + C.transpose()) / 5

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G * (Gp @ Gp) + Gp * (Gp @ G) + Gp * (G @ Gp)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs * (Gp @ Gp) + Gp * (Gp @ Gs) + Gp * (Gs @ Gp)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J.multiply(Jd * Jd) + Jd.multiply(Jd * J) + Jd.multiply(J * Jd)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js.multiply(Jd * Jd) + Jd.multiply(Jd * Js) + Jd.multiply(Js * Jd)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.multiply(Jd * Gd) + J.multiply(Gd * Jd) + G.multiply(Jd * Jd)
        C = C + Jd.multiply(Jd * G) + Jd.multiply(Gd * J) + Gd.multiply(Jd * J)
        C = C + Jd.multiply(J * Gd) + Jd.multiply(G * Jd) + Gd.multiply(J * Jd)
        motif_adj_mat = (C + C.transpose()) / 5

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js.multiply(Jd * Gd) + Js.multiply(Gd * Jd) + Gs.multiply(Jd * Jd)
        C = C + Jd.multiply(Jd * Gs) + Jd.multiply(Gd * Js) + Gd.multiply(Jd * Js)
        C = C + Jd.multiply(Js * Gd) + Jd.multiply(Gs * Jd) + Gd.multiply(Js * Jd)
        motif_adj_mat = (C + C.transpose()) / 5

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G.multiply(Gp * Gp) + Gp.multiply(Gp * G) + Gp.multiply(G * Gp)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs.multiply(Gp * Gp) + Gp.multiply(Gp * Gs) + Gp.multiply(Gs * Gp)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M4(adj_mat, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M4.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      Jd = mcin._build_Jd(adj_mat)
      motif_adj_mat = Jd * (Jd @ Jd)

    if mam_weight_type == "mean":
      Jd = mcin._build_Jd(adj_mat)
      Gd = mcin._build_Gd(adj_mat)
      motif_adj_mat = (Jd * (Jd @ Gd) + Jd * (Gd @ Jd) + Gd * (Jd @ Jd)) / 6

    if mam_weight_type == "product":
      Gp = mcin._build_Gp(adj_mat)
      motif_adj_mat = Gp * (Gp @ Gp)

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      Jd = mcin._build_Jd(adj_mat)
      motif_adj_mat = Jd.multiply(Jd * Jd)

    if mam_weight_type == "mean":
      Jd = mcin._build_Jd(adj_mat)
      Gd = mcin._build_Gd(adj_mat)
      motif_adj_mat = (Jd.multiply(Jd * Gd) + Jd.multiply(Gd * Jd) + Gd.multiply(Jd * Jd)) / 6

    if mam_weight_type == "product":
      Gp = mcin._build_Gp(adj_mat)
      motif_adj_mat = Gp.multiply(Gp * Gp)

  return motif_adj_mat


def mam_M5(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M5.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        C = J * (J @ J) + J * (J @ J.transpose()) + J * (J.transpose() @ J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        C = Js * (Js @ Js) + Js * (Js @ Js.transpose())
        C += Js * (Js.transpose() @ Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (J @ G) + J * (G @ J) + G * (J @ J)
        C = C + J * (J @ G.transpose()) + J * (G @ J.transpose())
        C += G * (J @ J.transpose())
        C = C + J * (J.transpose() @ G) + J * (G.transpose() @ J)
        C += G * (J.transpose() @ J)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Js * (Js @ Gs) + Js * (Gs @ Js) + Gs * (Js @ Js)
        C = C + Js * (Js @ Gs.transpose())
        C += Js * (Gs @ Js.transpose()) + Gs * (Js @ Js.transpose())
        C = C + Js * (Js.transpose() @ Gs)
        C += Js * (Gs.transpose() @ Js) + Gs * (Js.transpose() @ Js)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        C = G * (G @ G) + G * (G @ G.transpose())
        C += G * (G.transpose() @ G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        C = Gs * (Gs @ Gs) + Gs * (Gs @ Gs.transpose())
        C += Gs * (Gs.transpose() @ Gs)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        C = J.multiply(J * J) + J.multiply(J * J.transpose()) + J.multiply(J.transpose() * J)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        C = Js.multiply(Js * Js) + Js.multiply(Js * Js.transpose())
        C += Js.multiply(Js.transpose() * Js)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.multiply(J * G) + J.multiply(G * J) + G.multiply(J * J)
        C = C + J.multiply(J * G.transpose()) + J.multiply(G * J.transpose())
        C += G.multiply(J * J.transpose())
        C = C + J.multiply(J.transpose() * G) + J.multiply(G.transpose() * J)
        C += G.multiply(J.transpose() * J)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Js.multiply(Js * Gs) + Js.multiply(Gs * Js) + Gs.multiply(Js * Js)
        C = C + Js.multiply(Js * Gs.transpose())
        C += Js.multiply(Gs * Js.transpose()) + Gs.multiply(Js * Js.transpose())
        C = C + Js.multiply(Js.transpose() * Gs)
        C += Js.multiply(Gs.transpose() * Js) + Gs.multiply(Js.transpose() * Js)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        C = G.multiply(G * G) + G.multiply(G * G.transpose())
        C += G.multiply(G.transpose() * G)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        C = Gs.multiply(Gs * Gs) + Gs.multiply(Gs * Gs.transpose())
        C += Gs.multiply(Gs.transpose() * Gs)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M6(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M6.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J * (J @ Jd)
        Cprime = Jd * (J.transpose() @ J)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js * (Js @ Jd)
        Cprime = Jd * (Js.transpose() @ Js)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (J @ Gd) + J * (G @ Jd) + G * (J @ Jd)
        Cprime = Jd * (J.transpose() @ G) + Jd * (G.transpose() @ J)
        Cprime += Gd * (J.transpose() @ J)
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js * (Js @ Gd) + Js * (Gs @ Jd) + Gs * (Js @ Jd)
        Cprime = Jd * (Js.transpose() @ Gs)
        Cprime += Jd * (Gs.transpose() @ Js) + Gd * (Js.transpose() @ Js)
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G * (G @ Gp)
        Cprime = Gp * (G.transpose() @ G)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs * (Gs @ Gp)
        Cprime = Gp * (Gs.transpose() @ Gs)
        motif_adj_mat = C + C.transpose() + Cprime

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J.multiply(J * Jd)
        Cprime = Jd.multiply(J.transpose() * J)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js.multiply(Js * Jd)
        Cprime = Jd.multiply(Js.transpose() * Js)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.multiply(J * Gd) + J.multiply(G * Jd) + G.multiply(J * Jd)
        Cprime = Jd.multiply(J.transpose() * G) + Jd.multiply(G.transpose() * J)
        Cprime += Gd.multiply(J.transpose() * J)
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js.multiply(Js * Gd) + Js.multiply(Gs * Jd) + Gs.multiply(Js * Jd)
        Cprime = Jd.multiply(Js.transpose() * Gs)
        Cprime += Jd.multiply(Gs.transpose() * Js) + Gd.multiply(Js.transpose() * Js)
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G.multiply(G * Gp)
        Cprime = Gp.multiply(G.transpose() * G)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs.multiply(Gs * Gp)
        Cprime = Gp.multiply(Gs.transpose() * Gs)
        motif_adj_mat = C + C.transpose() + Cprime

  return motif_adj_mat


def mam_M7(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M7.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J * (Jd @ J)
        Cprime = Jd * (J @ J.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js * (Jd @ Js)
        Cprime = Jd * (Js @ Js.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (Jd @ G) + J * (Gd @ J) + G * (Jd @ J)
        Cprime = Jd * (J @ G.transpose()) + Jd * (G @ J.transpose())
        Cprime += Gd * (J @ J.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js * (Jd @ Gs) + Js * (Gd @ Js) + Gs * (Jd @ Js)
        Cprime = Jd * (Js @ Gs.transpose()) + Jd * (Gs @ Js.transpose())
        Cprime += Gd * (Js @ Js.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G * (Gp @ G)
        Cprime = Gp * (G @ G.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs * (Gp @ Gs)
        Cprime = Gp * (Gs @ Gs.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = J.multiply(Jd * J)
        Cprime = Jd.multiply(J * J.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        C = Js.multiply(Jd * Js)
        Cprime = Jd.multiply(Js * Js.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J.multiply(Jd * G) + J.multiply(Gd * J) + G.multiply(Jd * J)
        Cprime = Jd.multiply(J * G.transpose()) + Jd.multiply(G * J.transpose())
        Cprime += Gd.multiply(J * J.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Js.multiply(Jd * Gs) + Js.multiply(Gd * Js) + Gs.multiply(Jd * Js)
        Cprime = Jd.multiply(Js * Gs.transpose()) + Jd.multiply(Gs * Js.transpose())
        Cprime += Gd.multiply(Js * Js.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = G.multiply(Gp * G)
        Cprime = Gp.multiply(G * G.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Gp = mcin._build_Gp(adj_mat)
        C = Gs.multiply(Gp * Gs)
        Cprime = Gp.multiply(Gs * Gs.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

  return motif_adj_mat


def mam_M8(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M8.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = J * (J @ Jn)
        Cprime = Jn * (J.transpose() @ J)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (Js @ J0)
        Cprime = J0 * (Js.transpose() @ Js)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (G @ Jn) + G * (J @ Jn)
        Cprime = Jn * (J.transpose() @ G) + Jn * (G.transpose() @ J)
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (Gs @ J0) + Gs * (Js @ J0)
        Cprime = J0 * (Js.transpose() @ Gs) + J0 * (Gs.transpose() @ Js)
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = G * (G @ Jn)
        Cprime = Jn * (G.transpose() @ G)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Gs * (Gs @ J0)
        Cprime = J0 * (Gs.transpose() @ Gs)
        motif_adj_mat = C + C.transpose() + Cprime

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(J, J) - J.multiply(J)
        Cprime = J.transpose() * J - Id.multiply(J.transpose() * J)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Js, Js) - Js.multiply(Js * Je)
        Cprime = Js.transpose() * Js - Je.multiply(Js.transpose() * Js)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(J, G) - J.multiply(G) + mcut._a_b_one(G, J) - G.multiply(J)
        Cprime = J.transpose() * G - Id.multiply(J.transpose() * G)
        Cprime = Cprime + G.transpose() * J - Id.multiply(G.transpose() * J)
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Js, Gs) - Js.multiply(Gs * Je)
        C = C + mcut._a_b_one(Gs, Js) - Gs.multiply(Js * Je)
        Cprime = Js.transpose() * Gs - Je.multiply(Js.transpose() * Gs)
        Cprime = Cprime + Gs.transpose() * Js - Je.multiply(Gs.transpose() * Js)
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(G, G) - G.multiply(G)
        Cprime = G.transpose() * G - Id.multiply(G.transpose() * G)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Gs, Gs) - Gs.multiply(Gs * Je)
        Cprime = Gs.transpose() * Gs - Je.multiply(Gs.transpose() * Gs)
        motif_adj_mat = C + C.transpose() + Cprime

  return motif_adj_mat


def mam_M9(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M9.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = J * (Jn @ J.transpose()) + Jn * (J @ J)
        C += J * (J.transpose() @ Jn)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (J0 @ Js.transpose()) + J0 * (Js @ Js)
        C += Js * (Js.transpose() @ J0)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (Jn @ G.transpose()) + G * (Jn @ J.transpose())
        C = C + Jn * (J @ G) + Jn * (G @ J)
        C = C + J * (G.transpose() @ Jn) + G * (J.transpose() @ Jn)
        motif_adj_mat = (C + C.transpose()) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (J0 @ Gs.transpose()) + Gs * (J0 @ Js.transpose())
        C = C + J0 * (Js @ Gs) + J0 * (Gs @ Js)
        C = C + Js * (Gs.transpose() @ J0) + Gs * (Js.transpose() @ J0)
        motif_adj_mat = (C + C.transpose()) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = G * (Jn @ G.transpose()) + Jn * (G @ G) + G * (G.transpose() @ Jn)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Gs * (J0 @ Gs.transpose()) + J0 * (Gs @ Gs)
        C += Gs * (Gs.transpose() @ J0)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(J, J.transpose()) - 2 * J.multiply(J.transpose()) + J * J
        C = C - Id.multiply(J * J) + mcut._a_b_one(J, J.transpose())
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Js, Js.transpose()) - Js.multiply(Je * Js.transpose())
        C = C + Js * Js - Je.multiply(Js * Js)
        C = C + mcut._a_b_one(Js, Js.transpose()) - Js.multiply(Js.transpose() * Je)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(J, G.transpose()) - 2 * J.multiply(G.transpose()) + J * G
        C = C + mcut._a_one_b(G, J.transpose()) - 2 * G.multiply(J.transpose()) + G * J
        C = C - Id.multiply(J * G) + mcut._a_b_one(J, G.transpose())
        C = C - Id.multiply(G * J) + mcut._a_b_one(G, J.transpose())
        motif_adj_mat = (C + C.transpose()) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Js, Gs.transpose()) - Js.multiply(Je * Gs.transpose())
        C = C + mcut._a_one_b(Gs, Js.transpose()) - Gs.multiply(Je * Js.transpose())
        C = C + Js * Gs - Je.multiply(Js * Gs)
        C = C + mcut._a_b_one(Js, Gs.transpose()) - Js.multiply(Gs.transpose() * Je)
        C = C + Gs * Js - Je.multiply(Gs * Js)
        C = C + mcut._a_b_one(Gs, Js.transpose()) - Gs.multiply(Js.transpose() * Je)
        motif_adj_mat = (C + C.transpose()) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(G, G.transpose()) - 2 * G.multiply(G.transpose()) + G * G
        C = C - Id.multiply(G * G) + mcut._a_b_one(G, G.transpose())
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Gs, Gs.transpose()) - Gs.multiply(Je * Gs.transpose())
        C = C + Gs * Gs - Je.multiply(Gs * Gs)
        C = C + mcut._a_b_one(Gs, Gs.transpose()) - Gs.multiply(Gs.transpose() * Je)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M10(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M10.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = J * (Jn @ J)
        Cprime = Jn * (J @ J.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (J0 @ Js)
        Cprime = J0 * (Js @ Js.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = J * (Jn @ G) + G * (Jn @ J)
        Cprime = Jn * (J @ G.transpose()) + Jn * (G @ J.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Js * (J0 @ Gs) + Gs * (J0 @ Js)
        Cprime = J0 * (Js @ Gs.transpose()) + J0 * (Gs @ Js.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = G * (Jn @ G)
        Cprime = Jn * (G @ G.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Gs * (J0 @ Gs)
        Cprime = J0 * (Gs @ Gs.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(J, J) - J.multiply(J)
        Cprime = J * J.transpose() - Id.multiply(J * J.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Js, Js) - Js.multiply(Je * Js)
        Cprime = Js * Js.transpose() - Je.multiply(Js * Js.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(J, G) - J.multiply(G) + mcut._a_one_b(G, J) - G.multiply(J)
        Cprime = J * G.transpose() - Id.multiply(J * G.transpose())
        Cprime = Cprime + G * J.transpose() - Id.multiply(G * J.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Js, Gs) - Js.multiply(Je * Gs)
        C = C + mcut._a_one_b(Gs, Js) - Gs.multiply(Je * Js)
        Cprime = Js * Gs.transpose() - Je.multiply(Js * Gs.transpose())
        Cprime = Cprime + Gs * Js.transpose() - Je.multiply(Gs * Js.transpose())
        motif_adj_mat = (C + C.transpose() + Cprime) / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_one_b(G, G) - G.multiply(G)
        Cprime = G * G.transpose() - Id.multiply(G * G.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_one_b(Gs, Gs) - Gs.multiply(Je * Gs)
        Cprime = Gs * Gs.transpose() - Je.multiply(Gs * Gs.transpose())
        motif_adj_mat = C + C.transpose() + Cprime

  return motif_adj_mat


def mam_M11(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements

  """
  Perform the motif adjacency matrix calculations for motif M11.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        J = mcin._build_J(adj_mat)
        C = Jd * (J @ Jn) + Jn * (Jd @ J) + J * (Jd @ Jn)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Js = mcin._build_Js(adj_mat)
        C = Jd * (Js @ J0) + J0 * (Jd @ Js) + Js * (Jd @ J0)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = Jd * (G @ Jn) + Gd * (J @ Jn)
        C = C + Jn * (Jd @ G) + Jn * (Gd @ J)
        C = C + J * (Gd @ Jn) + G * (Jd @ Jn)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Jd * (Gs @ J0) + Gd * (Js @ J0)
        C = C + J0 * (Jd @ Gs) + J0 * (Gd @ Js)
        C = C + Js * (Gd @ J0) + Gs * (Jd @ J0)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = Gp * (G @ Jn) + Jn * (Gp @ G) + G * (Gp @ Jn)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Gp * (Gs @ J0) + J0 * (Gp @ Gs) + Gs * (Gp @ J0)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        J = mcin._build_J(adj_mat)
        C = mcut._a_b_one(Jd, J) - Jd.multiply(J)
        C = C + Jd * J - Id.multiply(Jd * J)
        C = C + mcut._a_b_one(J, Jd) - J.multiply(Jd)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Js = mcin._build_Js(adj_mat)
        C = mcut._a_b_one(Jd, Js) - Jd.multiply(Js * Je)
        C = C + Jd * Js - Je.multiply(Jd * Js)
        C = C + mcut._a_b_one(Js, Jd) - Js.multiply(Jd * Je)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = mcut._a_b_one(Jd, G) - Jd.multiply(G) + mcut._a_b_one(Gd, J) - Gd.multiply(J)
        C = C + Jd * G - Id.multiply(Jd * G) + Gd * J - Id.multiply(Gd * J)
        C = C + mcut._a_b_one(J, Gd) - J.multiply(Gd) + mcut._a_b_one(G, Jd) - G.multiply(Jd)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = mcut._a_b_one(Jd, Gs) - Jd.multiply(Gs * Je)
        C = C + mcut._a_b_one(Gd, Js) - Gd.multiply(Js * Je)
        C = C + Jd * Gs - Je.multiply(Jd * Gs) + Gd * Js - Je.multiply(Gd * Js)
        C = C + mcut._a_b_one(Js, Gd) - Js.multiply(Gd * Je)
        C = C + mcut._a_b_one(Gs, Jd) - Gs.multiply(Jd * Je)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Id = mcin._build_Id(adj_mat)
        G = mcin._build_G(adj_mat)
        C = mcut._a_b_one(Gp, G) - Gp.multiply(G)
        C = C + Gp * G - Id.multiply(Gp * G)
        C = C + mcut._a_b_one(G, Gp) - G.multiply(Gp)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = mcut._a_b_one(Gp, Gs) - Gp.multiply(Gs * Je)
        C = C + Gp * Gs - Je.multiply(Gp * Gs)
        C = C + mcut._a_b_one(Gs, Gp) - Gs.multiply(Gp * Je)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M12(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements

  """
  Perform the motif adjacency matrix calculations for motif M12.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        J = mcin._build_J(adj_mat)
        C = Jd * (Jn @ J) + Jn * (J @ Jd) + J * (Jn @ Jd)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Js = mcin._build_Js(adj_mat)
        C = Jd * (J0 @ Js) + J0 * (Js @ Jd) + Js * (J0 @ Jd)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        J = mcin._build_J(adj_mat)
        C = Jd * (Jn @ G) + Gd * (Jn @ J)
        C = C + Jn * (J @ Gd) + Jn * (G @ Jd)
        C = C + J * (Jn @ Gd) + G * (Jn @ Jd)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Jd * (J0 @ Gs) + Gd * (J0 @ Js)
        C = C + J0 * (Js @ Gd) + J0 * (Gs @ Jd)
        C = C + Js * (J0 @ Gd) + Gs * (J0 @ Jd)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = Gp * (Jn @ G) + Jn * (G @ Gp) + G * (Jn @ Gp)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = Gp * (J0 @ Gs) + J0 * (Gs @ Gp) + Gs * (J0 @ Gp)
        motif_adj_mat = C + C.transpose()

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        J = mcin._build_J(adj_mat)
        C = mcut._a_one_b(Jd, J) - Jd.multiply(J)
        C = C + J * Jd - Id.multiply(J * Jd)
        C = C + mcut._a_one_b(J, Jd) - J.multiply(Jd)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Js = mcin._build_Js(adj_mat)
        C = mcut._a_one_b(Jd, Js) - Jd.multiply(Je * Js)
        C = C + Js * Jd - Je.multiply(Js * Jd)
        C = C + mcut._a_one_b(Js, Jd) - Js.multiply(Je * Jd)
        motif_adj_mat = C + C.transpose()

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        C = mcut._a_one_b(Jd, G) - Jd.multiply(G) + mcut._a_one_b(Gd, J) - Gd.multiply(J)
        C = C + J * Gd - Id.multiply(J * Gd) + G * Jd - Id.multiply(G * Jd)
        C = C + mcut._a_one_b(J, Gd) - J.multiply(Gd) + mcut._a_one_b(G, Jd) - G.multiply(Jd)
        motif_adj_mat = (C + C.transpose()) / 3

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = mcut._a_one_b(Jd, Gs) - Jd.multiply(Je * Gs)
        C = C + mcut._a_one_b(Gd, Js) - Gd.multiply(Je * Js)
        C = C + Js * Gd - Je.multiply(Js * Gd) + Gs * Jd - Je.multiply(Gs * Jd)
        C = C + mcut._a_one_b(Js, Gd) - Js.multiply(Je * Gd)
        C = C + mcut._a_one_b(Gs, Jd) - Gs.multiply(Je * Jd)
        motif_adj_mat = (C + C.transpose()) / 3

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Id = mcin._build_Id(adj_mat)
        G = mcin._build_G(adj_mat)
        C = mcut._a_one_b(Gp, G) - Gp.multiply(G)
        C = C + G * Gp - Id.multiply(G * Gp)
        C = C + mcut._a_one_b(G, Gp) - G.multiply(Gp)
        motif_adj_mat = C + C.transpose()

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        Je = mcin._build_Je(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        C = mcut._a_one_b(Gp, Gs) - Gp.multiply(Je * Gs)
        C = C + Gs * Gp - Je.multiply(Gs * Gp)
        C = C + mcut._a_one_b(Gs, Gp) - Gs.multiply(Je * Gp)
        motif_adj_mat = C + C.transpose()

  return motif_adj_mat


def mam_M13(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif M13.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jd * (Jd @ Jn)
        Cprime = Jn * (Jd @ Jd)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Jd * (Jd @ J0)
        Cprime = J0 * (Jd @ Jd)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Jd * (Gd @ Jn) + Gd * (Jd @ Jn) + Jn * (Jd @ Gd)
        motif_adj_mat = (C + C.transpose()) / 4

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        C = Jd * (Gd @ J0) + Gd * (Jd @ J0) + J0 * (Jd @ Gd)
        motif_adj_mat = (C + C.transpose()) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Gp * (Gp @ Jn)
        Cprime = Jn * (Gp @ Gp)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = Gp * (Gp @ J0)
        Cprime = J0 * (Gp @ Gp)
        motif_adj_mat = C + C.transpose() + Cprime

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(Jd, Jd) - Jd.multiply(Jd)
        Cprime = Jd * Jd - Id.multiply(Jd * Jd)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Jd, Jd) - Jd.multiply(Jd * Je)
        Cprime = Jd * Jd - Je.multiply(Jd * Jd)
        motif_adj_mat = C + C.transpose() + Cprime

    if mam_weight_type == "mean":
      if motif_type == "func":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(Jd, Gd) - Jd.multiply(Gd) + mcut._a_b_one(Gd, Jd) - Gd.multiply(Jd)
        C = C + Jd * Gd - Id.multiply(Jd * Gd)
        motif_adj_mat = (C + C.transpose()) / 4

      if motif_type == "struc":
        Jd = mcin._build_Jd(adj_mat)
        Gd = mcin._build_Gd(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Jd, Gd) - Jd.multiply(Gd * Je)
        C = C + mcut._a_b_one(Gd, Jd) - Gd.multiply(Jd * Je)
        C = C + Jd * Gd - Je.multiply(Jd * Gd)
        motif_adj_mat = (C + C.transpose()) / 4

    if mam_weight_type == "product":
      if motif_type == "func":
        Gp = mcin._build_Gp(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = mcut._a_b_one(Gp, Gp) - Gp.multiply(Gp)
        Cprime = Gp * Gp - Id.multiply(Gp * Gp)
        motif_adj_mat = C + C.transpose() + Cprime

      if motif_type == "struc":
        Gp = mcin._build_Gp(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = mcut._a_b_one(Gp, Gp) - Gp.multiply(Gp * Je)
        Cprime = Gp * Gp - Je.multiply(Gp * Gp)
        motif_adj_mat = C + C.transpose() + Cprime

  return motif_adj_mat


def mam_Mcoll(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif Mcoll.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jn * (J @ J.transpose())
        motif_adj_mat = C

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Js @ Js.transpose())
        motif_adj_mat = C

    if mam_weight_type == "mean":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jn * (J @ G.transpose()) + Jn * (G @ J.transpose())
        motif_adj_mat = C / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Js @ Gs.transpose()) + J0 * (Gs @ Js.transpose())
        motif_adj_mat = C / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jn * (G @ G.transpose())
        motif_adj_mat = C

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Gs @ Gs.transpose())
        motif_adj_mat = C

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = J * J.transpose() - Id.multiply(J * J.transpose())
        motif_adj_mat = C

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Js * Js.transpose() - Je.multiply(Js * Js.transpose())
        motif_adj_mat = C

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = J * G.transpose() - Id.multiply(J * G.transpose())
        C = C + G * J.transpose() - Id.multiply(G * J.transpose())
        motif_adj_mat = C / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Js * Gs.transpose() - Je.multiply(Js * Gs.transpose())
        C = C + Gs * Js.transpose() - Je.multiply(Gs * Js.transpose())
        motif_adj_mat = C / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = G * G.transpose() - Id.multiply(G * G.transpose())
        motif_adj_mat = C

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Gs * Gs.transpose() - Je.multiply(Gs * Gs.transpose())
        motif_adj_mat = C

  return motif_adj_mat


def mam_Mexpa(adj_mat, motif_type, mam_weight_type, mam_method):

  """
  Perform the motif adjacency matrix calculations for motif Mexpa.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix from which to build the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to build.
  mam_weight_type : str
    The weighting scheme to use. One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    Which formulation to use. One of `"dense"` or `"sparse"`.

  Returns
  -------
  sparse matrix
    A motif adjacency matrix.
  """

  if mam_method == "dense":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jn * (J.transpose() @ J)
        motif_adj_mat = C

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Js.transpose() @ Js)
        motif_adj_mat = C

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        G = mcin._build_G(adj_mat)
        C = Jn * (J.transpose() @ G) + Jn * (G.transpose() @ J)
        motif_adj_mat = C / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Js.transpose() @ Gs) + J0 * (Gs.transpose() @ Js)
        motif_adj_mat = C / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Jn = mcin._build_Jn(adj_mat)
        C = Jn * (G.transpose() @ G)
        motif_adj_mat = C

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        J0 = mcin._build_J0(adj_mat)
        C = J0 * (Gs.transpose() @ Gs)
        motif_adj_mat = C

  if mam_method == "sparse":
    if mam_weight_type == "unweighted":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = J.transpose() * J - Id.multiply(J.transpose() * J)
        motif_adj_mat = C

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Js.transpose() * Js - Je.multiply(Js.transpose() * Js)
        motif_adj_mat = C

    if mam_weight_type == "mean":
      if motif_type == "func":
        J = mcin._build_J(adj_mat)
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = J.transpose() * G - Id.multiply(J.transpose() * G)
        C = C + G.transpose() * J - Id.multiply(G.transpose() * J)
        motif_adj_mat = C / 2

      if motif_type == "struc":
        Js = mcin._build_Js(adj_mat)
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Js.transpose() * Gs - Je.multiply(Js.transpose() * Gs)
        C = C + Gs.transpose() * Js - Je.multiply(Gs.transpose() * Js)
        motif_adj_mat = C / 2

    if mam_weight_type == "product":
      if motif_type == "func":
        G = mcin._build_G(adj_mat)
        Id = mcin._build_Id(adj_mat)
        C = G.transpose() * G - Id.multiply(G.transpose() * G)
        motif_adj_mat = C

      if motif_type == "struc":
        Gs = mcin._build_Gs(adj_mat)
        Je = mcin._build_Je(adj_mat)
        C = Gs.transpose() * Gs - Je.multiply(Gs.transpose() * Gs)
        motif_adj_mat = C

  return motif_adj_mat
