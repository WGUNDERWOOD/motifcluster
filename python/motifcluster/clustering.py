"""
Functions for spectral clustering
are in `motifcluster.clustering`.
"""

from sklearn.cluster import KMeans

from motifcluster import spectral as mcsp

def cluster_spectrum(spectrum, num_clusts):

  """
  Get cluster assignments from spectrum using k-means++.

  Get a list of cluster assignments from a spectrum,
  using k-means++ and `num_clusts` clusters.

  Parameters
  ----------
  spectrum : dict
    A dictionary containing `"vects"`:
    the matrix of eigenvectors to pass to k-means++.
  num_clusts : int
    The number of clusters to find.

  Returns
  -------
  cluster_assigns : list of int
    A list of integers from `1` to `num_clusts`,
    representing cluster assignments.
  """

  vects = spectrum["vects"][:, 1:]
  kmeans_plus_plus = KMeans(n_clusters=num_clusts).fit(vects)
  cluster_assigns = kmeans_plus_plus.labels_

  return cluster_assigns


def run_motif_clustering(adj_mat, motif_name,
                         motif_type="struc",
                         mam_weight_type="unweighted",
                         mam_method="sparse",
                         num_eigs=2,
                         type_lap="comb",
                         num_clusts=2,
                         restrict=True,
                         gr_method="sparse"):

  """
  Run motif-based clustering.

  Run motif-based clustering on the adjacency matrix of a
  (weighted directed) network,
  using a specified motif, motif type, weighting scheme,
  embedding dimension, number of clusters and Laplacian type.
  Optionally restrict to the largest connected component
  before clustering.

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
  num_clusts : int
    The number of clusters to find.
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
  vects : matrix
    A matrix containing the eigenvectors associated with the Laplace
    embedding of the (restricted) motif adjacency matrix.
  clusts :
    A vector containing integers representing the
    cluster assignment of each vertex in the (restricted) graph.

  Examples
  --------
  >>> adj_mat = np.array(range(1, 10)).reshape((3, 3))
  >>> run_motif_clustering(adj_mat, "M1")
  """

  assert motif_type in ["struc", "func"]
  assert mam_weight_type in ["unweighted", "mean", "product"]
  assert mam_method in ["sparse", "dense"]
  assert isinstance(restrict, bool)
  assert type_lap in ["comb", "rw"]

  spectrum = mcsp.run_motif_embedding(adj_mat, motif_name, motif_type, mam_weight_type,
                                      mam_method, num_eigs, type_lap, restrict, gr_method)

  cluster_assigns = cluster_spectrum(spectrum, num_clusts)

  ans = {
    "adj_mat": adj_mat,
    "motif_adj_mat": spectrum["motif_adj_mat"],
    "comps": spectrum["comps"],
    "adj_mat_comps": spectrum["adj_mat_comps"],
    "motif_adj_mat_comps": spectrum["motif_adj_mat_comps"],
    "vals": spectrum["vals"],
    "vects": spectrum["vects"],
    "clusts": cluster_assigns
  }

  return ans
