using TreeKnit
using TreeTools
using ArgParse
using StatsBase 
using MTKTools


function write_output_to_file(outfolder, filename, full_data; k_range=[0,1, 2])
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    open(outfolder*"/"*filename, "w") do io
        write(io, "k\tlCI\tuCI\tmean\tmedian\n")
        for k in k_range
            data = full_data[k]
            write(io, string(k) *"\t"*string(data["lCI"])*"\t"*string(data["uCI"])*"\t"*string(data["mean"])*"\t"*string(data["median"])*"\n")
        end
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
        "--sim"
            help = "number of simulations"
            arg_type = Int
            default = 10
        "--simtype"
            help = "simulate under flu or kingman model"
            arg_type = String
            default = "flu"
        "--res"
            help = "rate of resolution"
            arg_type = Float64
            default = 0.3
    end

    return parse_args(s)
end

function add_to_dict_vec(dict_, key, value)
    if !haskey(dict_, key)
        dict_[key] = [value]
    else
        append!(dict_[key], value)
    end
end

function get_percentage_ambiguous_splits(t1, t2, MCCs)
    resolvable_splits_strict = TreeKnit.new_splits(MCCs, t1, t2; strict=true)
    resolvable_splits_liberal = TreeKnit.new_splits(MCCs, t1, t2; strict=false)
    println("splits")
    println(resolvable_splits_strict[1].splits)
    println(resolvable_splits_liberal[1].splits)
    return (length(resolvable_splits_liberal[1].splits) - length(resolvable_splits_strict[1].splits)), length(resolvable_splits_strict[1].splits)
end

function run_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; 
                        simtype = :flu, res=0.3)
    r = 10^rec_rate
    ambiguous_splits = Dict{Int, Vector{Float64}}()   
    
    for i in 1:no_sim
        c = get_c(res, rec_rate; n=no_lineages, simtype)
        true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)
        rMCCs = MTKTools.get_real_MCCs(8, arg)
        rMCCs = TreeKnit.MCC_set(8, [t.label for t in true_trees], rMCCs)

        unresolved_trees = [copy(t) for t in true_trees]
        unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

        rand_MCC = sample(1:8, 2, replace = false)
        MCCs = get(rMCCs, rand_MCC...)
        t1 = unresolved_trees[rand_MCC[1]]
        t2 = unresolved_trees[rand_MCC[2]]
        average_size_MCC = sum([length(m) for m in MCCs])/length(MCCs)

        ambig, non_ambig_splits = get_percentage_ambiguous_splits(t1, t2, MCCs)
        add_to_dict_vec(ambiguous_splits, 0, average_size_MCC)

        add_to_dict_vec(ambiguous_splits, 1, ambig)
        add_to_dict_vec(ambiguous_splits, 2, non_ambig_splits)
    end

    ambiguous_splits_vector = Dict()
    for i in [0,1, 2]
        cs_ = sort(ambiguous_splits[i])
        CI_ = [cs_[Int(round(no_sim*0.05))+1], cs_[Int(round(no_sim*0.95))]]
        correct_data = Dict("lCI" =>CI_[1], "mean" => sum(ambiguous_splits[i])/no_sim, "uCI" =>CI_[2], "median" => cs_[Int(round(no_sim*0.5))] )

        ambiguous_splits_vector[i] = correct_data
    end  
    return  ambiguous_splits_vector
end

function main()
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    num_sim = parsed_args["sim"]
    simt = parsed_args["simtype"]
    res = parsed_args["res"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    results = run_simulations(num_sim, n, r; simtype, res)
    filename = "results_"*string(r)*".txt"
    write_output_to_file(o, filename, results; k_range=[0,1, 2])
end

main()