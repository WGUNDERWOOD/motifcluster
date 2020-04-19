cluster_spectrum <- function(spectrum, num_clusts) {

  vects = spectrum$vects[, -1]
  kmeans_plus_plus <- LICORS::kmeanspp(vects, num_clusts)
  cluster_assigns <- kmeans_plus_plus$cluster

  return(cluster_assigns)
}
