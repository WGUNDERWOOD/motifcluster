using Graphs
using MotifCluster

function performance_trial(ns::Vector{Int}, k::Int, motifs::Vector{String},
        nreps::Int, graph_type::String)

    results = Dict{String, Any}[]
    for n in sort(ns, rev=true)
        for motif in motifs
            for rep in 1:nreps

                # graph parameters
                block_sizes = [n]
                connection_matrix = reshape([k/n], (1, 1))

                # sample graph
                if graph_type == "erdos_renyi"
                    adj_mat = sample_dsbm(block_sizes, connection_matrix, nothing, "unweighted")
                elseif graph_type == "barabasi_albert"
                    sample_graph = barabasi_albert(n, k)
                    adj_mat = adjacency_matrix(sample_graph)
                end

                # time mam construction
                t0 = time()
                build_motif_adjacency_matrix(adj_mat, motif, "func", "mean")
                t1 = time()
                dt = t1 - t0

                # add results
                result = Dict("n" => n, "k" => k, "motif" => motif, "rep" => rep, "time" => dt)
                push!(results, result)
                println("n = ", n, ", k = ", k, ", ", motif, ", rep = ", rep, ", graph type = ",
                        graph_type, ", time = ", round(dt, digits=3))

            end
        end
    end

    # save results
    filename = "results/julia_k" * string(k) * "_" * graph_type * ".csv"
    open(filename, "w") do file
        write(file, "n,k,motif,time,rep\n")
        for result in results
            write(file, string(result["n"]) * ",")
            write(file, string(k) * ",")
            write(file, result["motif"] * ",")
            write(file, string(result["time"]) * ",")
            write(file, string(result["rep"]) * "\n")
        end
    end

    return nothing
end

motifs = ["M1","M8","M11"]
nreps = 10

ns= [101, 200, 500, 1000]
performance_trial(ns, 100, motifs, nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, nreps, "erdos_renyi")

ns= [101, 200, 500, 1000]
performance_trial(ns, 100, motifs, nreps, "barabasi_albert")
performance_trial(ns, 10, motifs, nreps, "barabasi_albert")
performance_trial(ns, 100, motifs, nreps, "erdos_renyi")
performance_trial(ns, 10, motifs, nreps, "erdos_renyi")
