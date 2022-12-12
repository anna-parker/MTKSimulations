using TreeAlgs.Evolve
using TreeTools
using CSV
using TreeKnit
using ArgParse
using Random
using Combinatorics
using Distributions
using MTKTools

function write_real_tree!(tree::Tree, μ::Number, tree_name::String, outfolder::AbstractString)
    #normalize tree by mutation rate
    for node in POT(tree)
        node.tau= node.tau * μ
    end

    write_newick(outfolder *"/"*tree_name *".nwk", tree)
end

"""
    simulate(n::Int, r::Number, outfolder::AbstractString; flu = true, debug=true)
Simulate a pair of ARG trees using sample size (`n`), and recombination rate ratio to coalescence (`r`),
mutate sequences along simulated trees (if `flu` use sequence lengths of HA and NA segments, else 1000).
"""
function simulate(n::Int, rec_rate::Number, outfolder::AbstractString; remove=false, rand_error=true, infer=false, res_rate=0.3, simtype = :flu, strict=true, s=0.0)
    ##set sequence length length and mutation rate
    r = 10^rec_rate
    L = 1000
    N = 10_000
    c = MTKTools.get_c(res_rate, rec_rate; n=n, simtype)
    μ = 4/(N*c*L*3)

    ##simulate trees and evolve segments
    trees, arg = MTKTools.get_trees(8, n; c, ρ = r, simtype, s)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    for i in 1:8
        trees[i].label = string(i)
        evolve!(trees[i], L, μ);
        write_seq2fasta(trees[i], string(i), outfolder, only_terminals=true, remove_0_mutations=false, write_date=false);
        write_real_tree!(trees[i], μ, "true_segment_"*string(i), outfolder)
    end

    ## create metadata files with dates Int(
    node_labels = String[]
    for node in POTleaves(trees[1])
        push!(node_labels, node.label)
    end
    start_dates = fill("2000.01" , length(node_labels))
    CSV.write(outfolder *"/metadata.csv", (name=node_labels, date=start_dates))

    delete_list = Vector{String}[]
    for tree in trees
        l = []
        for node in internals(tree)
            if !node.isroot
                if !node.data.dat["evolved"]
                    push!(l, node.label)
                end
            end
        end
        append!(delete_list, [l])
    end
    print("Average res rate: "*string(((n-1-mean([length(l) for l in delete_list]))/(n-1)))*"\n")
    if remove
        ## remove internal nodes (except the root) with no mutations between them and their children
        #trees = remove_branches(trees; c=c, N = 10_000)
        for node_label in delete_list
            delete_node!(trees[i], node_label, ptau=true)
        end
    end
    if rand_error
        d = Normal(0, 0.05)
        td = truncated(d, -1, 1)
        ##write trees with relative error to nwk files
        for i in range(1,2)
            tree = trees[i]
            for node in POT(tree)
                if !node.isroot
                    node.tau = (1-rand(td))*node.tau 
                end
            end
        end
    end
    rMCCs = MTKTools.get_real_MCCs(8, arg)
    rMCCs = TreeKnit.MCC_set(8, [string(i) for i in 1:8], rMCCs)
    for no_trees in 2:8
        rand_order = sample(1:8, no_trees, replace = false)
        i_trees = [copy(t) for t in trees[rand_order]]
        if infer == true
            i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(;nMCMC=250, consistent = false, parallel=true, strict))
        else
            i_MCCs = Vector{Vector{String}}[]
            for pair in combinations(rand_order,2)
                append!(i_MCCs, [TreeKnit.get(rMCCs, pair...)])
            end
            i_MCCs = TreeKnit.MCC_set(no_trees, [string(r) for r in rand_order], i_MCCs)
        end
        if !isdir(outfolder * "/"* string(no_trees))
            mkdir(outfolder * "/"* string(no_trees))
        end
        write_mccs(outfolder * "/"* string(no_trees)* "/" * "MCCs.json", i_MCCs)
    end

    ## write modified tree to nwk files
    mkpath(outfolder)
    for i in 1:8
        write_newick(outfolder *"/"*string(i)*".nwk", trees[i])
    end

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--n"
            help = "sample size at time 0 (Int)"
            arg_type = Int
            default = 50
        "--rec"
            help = "recombination rate ratio to coalescence (Float)"
            arg_type = Float64
            default = 0.0
        "--o"
            help = "output directory"
            arg_type = String
            default = "results"
        "--i"
            help = "if inferred or real MCCs"
            arg_type = String
            default = "infer"
        "--simtype"
            help = "which simtype to use"
            arg_type = String
            default = "flu"
        "--res"
            help = "resolution rate"
            arg_type = Float64
            default = 0.3
        "--strict"
            help = "strict or liberal resolve"
            arg_type = Bool
            default = true
        "--s"
            help= "rate to add nodes backwards in time"
            arg_type = Float64
            default = 0.0
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    i = parsed_args["i"]
    simt = parsed_args["simtype"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    res = parsed_args["res"]
    strict = parsed_args["strict"]
    s = parsed_args["s"]
    # n = 10
    # i = "real"
    # r = -4.0
    # o = "DivergenceTimeEstimations/sim_test/sim_1_1"
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r")
    if i=="infer"
        simulate(n, r, o; res_rate=res, simtype, strict, s)
    else
        simulate(n, r, o; infer=false, res_rate=res, simtype, strict, s)
    end
end

main()