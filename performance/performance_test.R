library(devtools)
library(igraph)
load_all("../R/")



# functions
#######################################################################

performance_trial = function(ns, k, motifs, method, nreps, graph_type){

  results = data.frame(
    n = integer(),
    k = double(),
    motif = character(),
    method = character(),
    time = double()
  )

  for(n in sort(ns, decreasing = TRUE)){
    for(motif in motifs){
      for(rep in 1:nreps){

        # graph parameters
        block_sizes = c(n)
        connection_matrix = matrix(k / n)

        # sample graph
        if(graph_type == "erdos_renyi"){
          adj_mat = sample_dsbm(block_sizes, connection_matrix)
        }
        else if(graph_type == "barabasi_albert"){
          sample_graph = sample_pa(
            n, k, directed = FALSE,
            start.graph = make_empty_graph(n = k, directed = FALSE)
          )
          adj_mat = drop0(as_adjacency_matrix(sample_graph))
        }

        if(method == "dense"){
          adj_mat = matrix(adj_mat, nrow = n)
        }

        # time mam construction
        t0 = Sys.time()
        build_motif_adjacency_matrix(adj_mat, motif, "func", "mean", method)
        t1 = Sys.time()
        dt = t1 - t0

        # add results
        results_row = data.frame("n" = n, "k" = k, "motif" = motif, "method" = method, "time" = dt, "rep" = rep)
        results = rbind(results, results_row)

        cat("n = ",n, ", k = ",k, ", ", motif, ", method = ", method,
              ", rep = ", rep, ", graph type = ", graph_type, ", time = ", round(dt, 3), "\n", sep = "")
        }
      }
    }

  # save results
  csv_name = paste("results/r_k", k, "_", method, "_", graph_type, ".csv", sep = "")
  write.csv(results, csv_name, quote = FALSE, row.names = FALSE)
}



# script
#######################################################################

motifs = c('M1','M8','M11')
nreps = 2

ns= c(101, 200, 500, 1000)
performance_trial(ns, 100, motifs, "dense", nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, "dense", nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, "dense", nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, "dense", nreps, "erdos_renyi")

ns= c(101, 200, 500, 1000, 2000, 5000)
performance_trial(ns, 100, motifs, "sparse", nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, "sparse", nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, "sparse", nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, "sparse", nreps, "erdos_renyi")
