push!(LOAD_PATH, "../GtensorSC")
include("../GtensorSC/util.jl")

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--alpha", "-a"
            help = "alpha parameter"
            arg_type = Float64
            default = 0.95
        "--phi", "-p"
            help = "phi conductance threshold"
            arg_type = Float64
            default = 0.8
        "--min-cut"
            help = "minimum cluster size"
            arg_type = Int
            default = 4
        "--max-cut"
            help = "maximum cluster size"
            arg_type = Int
            default = 100
        "infile"
            help = "file containing the tensor"
            required = true
        "outfile"
            help = "file to write the clusters"
            required = true
    end

    return parse_args(s)
end

# Read tensor
args = parse_commandline()
filename = args["infile"]
# filename = "GtensorSC/data/openFlight/tensor_openflight_small.txt"
# filename = "data/net_merge_final.txt"
P = read_tensor(filename);

# Set parameters, run algorithm
para = algPara(args["alpha"], args["min-cut"], args["max-cut"], args["phi"])
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

open(args["outfile"], "w") do f
        for (cluster, members) in clusters
            write(f, string("#cluster ", cluster, "\n"))
            write(f, string(join(members, " "), "\n"))
        end
    end
