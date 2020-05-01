"""
Functions relating to spectral methods
are in `motifcluster.spectral`.
"""

import numpy as np
from scipy import sparse
from scipy import linalg

from motifcluster import utils as mcut
from motifcluster import motifadjacency as mcmo

def _get_first_eigs(some_mat, num_eigs):

  """
  Compute first few eigenvalues and eigenvectors of a matrix.

  Compute the first few eigenvalues (by magnitude) and
  associated eigenvectors of a matrix.

  Parameters
  ----------
  some_mat : matrix
    Matrix for which eigenvalues and eigenvectors
    are to be calculated.
  num_eigs : int
    Number of eigenvalues and eigenvectors to calculate.

  Returns
  -------
  vals : list
    A length-`num_eigs` list of the first few eigenvalues.
  vects : matrix
    A `some_mat.shape[0]` by `num_eigs` matrix
    of the associated eigenvectors.
  """

  # check args
  assert num_eigs == np.floor(num_eigs)
  assert num_eigs >= 1

  # get spectrum for many eigs
  if num_eigs >= some_mat.shape[0] - 1:

    # make sure matrix is dense
    if sparse.issparse(some_mat):
      some_mat = some_mat.toarray()

    vals, vects = linalg.eig(some_mat)
    ordering = np.argsort(vals.real)
    vals = vals.real[ordering]
    vects = vects.real[:, ordering]
    vals = vals[0:num_eigs]
    vects = vects[:, 0:num_eigs]

  # get spectrum for few eigs
  else:

    # make sure matrix is sparse
    if not sparse.issparse(some_mat):
      some_mat = sparse.csr_matrix(some_mat, dtype="f")

    vals, vects = sparse.linalg.eigs(some_mat, k=num_eigs, which="SM")
    ordering = np.argsort(vals.real)
    vals = vals.real[ordering]
    vects = vects.real[:, ordering]

  # return a list
  spectrum = {
    "vects": vects,
    "vals": vals
  }

  return spectrum


def build_laplacian(adj_mat, type_lap="rw"):

  """
  Build a Laplacian matrix.

  Build a Laplacian matrix (combinatorial Laplacian or random-walk Laplacian)
  from a symmetric (weighted) graph adjacency matrix.

  Parameters
  ----------
  adj_mat : matrix
    Symmetric adjacency matrix from which to build the Laplacian.
  type_lap : str
    Type of Laplacian to build.
    One of `"comb"` (combinatorial) or `"rw"` (random-walk).

  Returns
  -------
  sparse matrix
    The specified Laplacian matrix.

  Examples
  --------
  >>> adj_mat = np.array(range(1, 10)).reshape((3, 3))
  >>> build_laplacian(adj_mat, "rw")
  """

  # check args
  assert type_lap in ["comb", "rw"]
  n = adj_mat.shape[0]

  if not sparse.issparse(adj_mat):
    adj_mat = sparse.csr_matrix(adj_mat)

  # initialize parameters
  degs_adj_mat = adj_mat.sum(axis=0).reshape(n)

  # combinatorial Laplacian
  if type_lap == "comb":
    degs_matrix = sparse.diags(degs_adj_mat, offsets=[0], shape=(n, n))
    L = degs_matrix - adj_mat

  # random-walk Laplacian
  elif type_lap == "rw":
    assert (degs_adj_mat > 0).all()
    inv_degs_matrix = sparse.diags(1 / degs_adj_mat, offsets=[0], shape=(n, n))
    L = sparse.identity(n) - inv_degs_matrix * adj_mat

  L = sparse.csr_matrix(L)

  return L


def run_laplace_embedding(adj_mat, num_eigs, type_lap="rw"):

  """
  Run Laplace embedding.

  Run Laplace embedding on a symmetric (weighted) adjacency matrix
  with a specified number of eigenvalues and eigenvectors.

  Parameters
  ----------
  adj_mat : matrix
    Symmetric adjacency matrix to be embedded.
  num_eigs : int
    Number of eigenvalues and eigenvectors for the embedding.
  type_lap : str
    Type of Laplacian for the embedding.
    One of `"comb"` (combinatorial) or `"rw"` (random-walk).

  Returns
  -------
  vals : list
    The length-`num_eigs` list
    of the first few eigenvalues of the Laplacian.
  vects : matrix
    An `adj_mat.shape[0]` by `num_eigs` matrix
    of the associated eigenvectors.

  Examples
  --------
  >>> adj_mat = np.array(range(1, 10)).reshape((3, 3))
  >>> run_laplace_embedding(adj_mat, 2, "rw")
  """

  # check args
  assert num_eigs == np.floor(num_eigs)
  assert num_eigs >= 1
  assert type_lap in ["comb", "rw"]

  # build and embed Laplacian
  laplacian = build_laplacian(adj_mat, type_lap)
  spectrum = _get_first_eigs(laplacian, num_eigs)

  return spectrum


def run_motif_embedding(adj_mat, motif_name,
                        motif_type="struc",
                        mam_weight_type="unweighted",
                        mam_method="sparse",
                        num_eigs=2,
                        type_lap="rw",
                        restrict=True,
                        gr_method="sparse"):
  """
  Run motif embedding.

  Calculate a motif adjacency matrix for a given motif and motif type,
  optionally restrict it to its largest connected component,
  and then run Laplace embedding with specified Laplacian type and
  number of eigenvalues and eigenvectors.

  Parameters
  ----------
  adj_mat : matrix
    Adjacency matrix to be embedded.
  motif_name : str
    Motif used for the motif adjacency matrix.
  motif_type : str
    Type of motif adjacency matrix to use.
    One of `"func"` or `"struc"`.
  mam_weight_type : str
    Weighting scheme for the motif adjacency matrix.
    One of `"unweighted"`, `"mean"` or `"product"`.
  mam_method : str
    The method to use for building the motif adjacency matrix.
    One of `"sparse"` or `"dense"`.
  num_eigs : int
    Number of eigenvalues and eigenvectors for the embedding.
  type_lap : str
    Type of Laplacian for the embedding.
    One of `"comb"` or `"rw"`.
  restrict : bool
    Whether or not to restrict the motif adjacency matrix
    to its largest connected component before embedding.
  gr_method : str
    Format to use for getting largest component.
    One of `"sparse"` or `"dense"`.

  Returns
  -------
  adj_mat : sparse matrix
    The original adjacency matrix.
  motif_adj_mat : sparse matrix
    The motif adjacency matrix.
  comps : list
    The indices of the largest connected component
    of the motif adjacency matrix
    (if restrict=True).
  adj_mat_comps : matrix
    The original adjacency matrix restricted
    to the largest connected component of the motif adjacency matrix
    (if restrict=True).
  motif_adj_mat_comps : matrix
    The motif adjacency matrix restricted
    to its largest connected component
    (if restrict=True).
  vals : list
    A length-`num_eigs` list containing the
    eigenvalues associated with the Laplace embedding
    of the (restricted) motif adjacency matrix.
  vects :
    A matrix
    containing the eigenvectors associated with the Laplace embedding
    of the (restricted) motif adjacency matrix.

  Examples
  --------
  adj_mat = np.array(range(1, 10)),reshape((3, 3))
  run_motif_embedding(adj_mat, "M1")
  """

  # check args
  assert motif_name in mcut.get_motif_names()
  assert motif_type in ["struc", "func"]
  assert num_eigs == np.floor(num_eigs)
  assert num_eigs >= 1
  assert mam_weight_type in ["unweighted", "mean", "product"]
  assert mam_method in ["sparse", "dense"]
  assert type_lap in ["comb", "rw"]
  assert isinstance(restrict, bool)
  assert gr_method in ["sparse", "dense"]

  if not sparse.issparse(adj_mat):
    adj_mat = sparse.csr_matrix(adj_mat)

  # build motif adjacency matrix
  motif_adj_mat = mcmo.build_motif_adjacency_matrix(adj_mat, motif_name,
                                                    motif_type, mam_weight_type,
                                                    mam_method)

  if restrict:

    # restrict to largest connected component
    comps = mcut.get_largest_component(motif_adj_mat, gr_method)

    adj_mat_comps = adj_mat[comps, :].tocsc()[:, comps]
    adj_mat_comps = sparse.csr_matrix(adj_mat_comps)

    motif_adj_mat_comps = motif_adj_mat[comps, :].tocsc()[:, comps]
    motif_adj_mat_comps = sparse.csr_matrix(motif_adj_mat_comps)

    # Laplace embedding restricted
    spect = run_laplace_embedding(motif_adj_mat_comps, num_eigs, type_lap)

  else:
    comps = None
    adj_mat_comps = None
    motif_adj_mat_comps = None

    # Laplace embedding unrestricted
    spect = run_laplace_embedding(motif_adj_mat, num_eigs, type_lap)


  # return list
  embedding = {
    "adj_mat": adj_mat,
    "motif_adj_mat": motif_adj_mat,
    "comps": comps,
    "adj_mat_comps": adj_mat_comps,
    "motif_adj_mat_comps": motif_adj_mat_comps,
    "vals": spect["vals"],
    "vects": spect["vects"],
  }

  return embedding
