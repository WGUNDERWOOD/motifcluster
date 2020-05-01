library(motifcluster)

t00 = Sys.time()
for (n in c(100, 200, 500, 1000, 2000)) {
  for (p in c(10/n, 100/n)){
    for (motif_name in c("Ms", "Md", "M1", "M9", "M11")){

      t0 = Sys.time()
      block_sizes = c(n)
      connection_matrix = matrix(p)
      weight_matrix = matrix(10)

      adj_mat = sample_dsbm(block_sizes, connection_matrix,
                            weight_matrix)

      build_motif_adjacency_matrix(adj_mat, motif_name)

      cat(n, p, motif_name, "\n")
      cat(round(Sys.time() - t0, 2), "\n\n")
    }
  }
}

cat("Total time: \n")
cat(round(Sys.time() - t00, 2)"\n")
