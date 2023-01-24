using TreeAlgs.Evolve
using TreeTools
using CSV
using TreeKnit
using ArgParse
using Random
using Combinatorics
using Distributions
using MTKTools
using JSON3

"""
create_accuracy_shared_branches_dict(tree, true_tree, MCC, rMCC)

Calculate if the branch corresponding to a node is a TP (true positive), FP, TN, and FN (false negative)
for shared branch inference. Label the true underlying tree `true_tree` (prior to branch removal)
using the real MCCs `MCC`, do the same with the infered, resolved tree `tree` (assuming TreeKnit received
an unresolved version of `true_tree` with removed branches and branches were then added to that tree using 
the inferred MCCs) and inferred MCCs `MCC`. For each branch that exists in both trees (use splitlist)
check if it has been correctly labeled as shared or not shared.
"""
function create_accuracy_shared_branches_dict(tree, true_tree, MCC, rMCC)
    accuracy_dict = Dict{String, String}()

    ##note the true tree should be fully resolved
    shared_branches_ = MTK.map_shared_branches(MCC, tree)
    real_shared_branches_ = MTK.map_shared_branches(rMCC, true_tree)
    true_splits = SplitList(true_tree)
    true_splits_dict = Dict()
    for n in nodes(true_tree)
        if !isleaf(n) && !isroot(n)
            split = true_splits.splitmap[n.label]
            true_splits_dict[split] = n
        end
    end
    s1 = SplitList(tree)
    for n in nodes(tree)
        node_in_true_tree = nothing
        if !isleaf(n) && !isroot(n)
            split = s1.splitmap[n.label]
            if split ∈ SplitList(true_tree)
                ##true split, branch exists in true tree
                node_in_true_tree = true_splits_dict[split]
            end
        end
        if isleaf(n)
            node_in_true_tree = true_tree.lnodes[n.label]
        end
        if !isnothing(node_in_true_tree)
            if shared_branches_[n.label]
                if real_shared_branches_[node_in_true_tree.label]
                    accuracy_dict[n.label] = "TP"
                else
                    accuracy_dict[n.label] = "FP"
                end
            else
                if real_shared_branches_[node_in_true_tree.label]
                    accuracy_dict[n.label] = "FN"
                else
                    accuracy_dict[n.label] = "TN"
                end
            end
        end
    end
    return accuracy_dict
end

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
function simulate(n::Int, rec_rate::Number, outfolder::AbstractString; remove=false, rand_error=true, infer=false, res_rate=0.3, simtype = :flu, strict=true, s=0.0, pair=false)
    ##set sequence length length and mutation rate
    r = 10^rec_rate
    L = 1000
    N = 10_000
    c = MTKTools.get_c(res_rate, rec_rate; n=n, simtype)
    μ = 1/(N*c*L)

    ##simulate trees and evolve segments
    trees, arg = MTKTools.get_trees(8, n; c, ρ = r, simtype, s)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    if pair
        k_vals, no_trees_ = 1:2, [2]
    else
        k_vals, no_trees_ =1:8, 2:8
    end
    for i in k_vals
        trees[i].label = string(i)
        evolve!(trees[i], L, (4/3)* μ);
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
    for tree in trees[k_vals]
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
        for (t, node_label) in (trees[k_vals], delete_list[k_vals])
            delete_node!(t, node_label, ptau=true)
        end
    end
    if rand_error
        d = Normal(0, 0.05)
        td = truncated(d, -1, 1)
        ##write trees with relative error to nwk files
        for t in trees[k_vals]
            for node in POT(t)
                if !node.isroot
                    node.tau = (1-rand(td))*node.tau 
                end
            end
        end
    end
    rMCCs = MTKTools.get_real_MCCs(8, arg)
    rMCCs = TreeKnit.MCC_set(8, [t.label for t in trees], rMCCs)
    for no_trees in no_trees_
        if pair
            rand_order = [1,2]
        else
            rand_order = sample(1:8, no_trees, replace = false)
        end
        i_trees = [copy(t) for t in trees[rand_order]]
        names = [t.label for t in i_trees]
        if infer == true
            i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(;nMCMC=250, parallel=true, rounds=1, resolve=false, strict))
        else
            i_MCCs = Vector{Vector{String}}[]
            for pair in combinations(rand_order,2)
                append!(i_MCCs, [TreeKnit.get(rMCCs, pair...)])
            end
            i_MCCs = TreeKnit.MCC_set(no_trees, [t.label for t in trees[rand_order]], i_MCCs)
        end
        if pair
            dict_ = create_accuracy_shared_branches_dict(i_trees[1], i_trees[1], TreeKnit.get(i_MCCs, names...), TreeKnit.get(rMCCs, names...))
            json_string = JSON3.write(dict_)
            open(outfolder * "/"* "shared_branches.json", "w") do f
                JSON3.pretty(f, json_string)
            end
        end
        if !isdir(outfolder * "/"* string(no_trees))
            mkdir(outfolder * "/"* string(no_trees))
        end
        write_mccs(outfolder * "/"* string(no_trees)* "/" * "MCCs.json", i_MCCs)
    end

    ## write modified tree to nwk files
    mkpath(outfolder)
    for i in k_vals
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
            default = -4.0
        "--o"
            help = "output directory"
            arg_type = String
            default = "results"
        "--i"
            help = "if inferred or real MCCs"
            arg_type = Bool
            default = true
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
        "--segments"
            help= "only look at a pair and write shared branches to json"
            arg_type = Int
            default = 8
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
    segments = parsed_args["segments"]
    if segments == 2
        pair = true
    else
        pair = false
    end
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r")
    simulate(n, r, o; infer=i, res_rate=res, simtype, strict, s, pair)
end

main()