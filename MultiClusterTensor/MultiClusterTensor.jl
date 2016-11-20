push!(LOAD_PATH, "../GtensorSC")
include("../GtensorSC/util.jl")

# Read tensor
filename = ARGS[1]
# filename = "GtensorSC/data/openFlight/tensor_openflight_small.txt"
# filename = "data/net_merge_final.txt"
P = read_tensor(filename);

# Set parameters, run algorithm
para = algPara(0.95, 5, 100, 0.8)
(r,h) = tensor_speclustering(P, para);

# Parse into clusters
# indVec is the clustering result indicating group number for each node
indVec = zeros(r.n)
traCount = trav_tree(r, indVec, 1)

# Save clusters as lines
clusters = Dict()
for node_i in 1:length(indVec)
    ind = round(Int, indVec[node_i])
    if haskey(clusters, ind)
        push!(clusters[ind], node_i)
    else
        clusters[ind] = [node_i]
    end
end

open(ARGS[2], "w") do f
        for (cluster, members) in clusters
            write(f, string("#cluster ", cluster, "\n"))
            write(f, string(join(members, " "), "\n"))
        end
    end
