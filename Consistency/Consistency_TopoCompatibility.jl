using MTKTools
using TreeTools
using TreeKnit
using TreeKnit.MTK
using TestRecombTools
using ArgParse
using StatsBase
using Clustering
using Dagger


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
        "--strict"
            help = "strict or liberal resolve"
            arg_type = Bool
            default = true
        "--rounds"
            help = "number of rounds"
            arg_type = Int
            default = 2
        "--krange"
            help = "if full 2:8 (true) or [2,4,8] (false)"
            arg_type = Bool
            default = false
        "--final-no-resolve"
            help = "default false"
            arg_type = Bool
            default = false
        "--pre-resolve"
            help = "default false"
            arg_type = Bool
            default = false
    end

    return parse_args(s)
end

function write_output(outfolder, data_standard, type, rec_rate; k_range=2:8)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    filename = "results_"*type*"_"*string(rec_rate)*".txt"
    write_output_to_file(outfolder*"/"*filename, data_standard; k_range)
end

function write_output_to_file(filename, full_data; k_range=2:8)
    open(filename, "w") do io
        write(io, "k\tlCI\tuCI\tmean\tmedian\n")
        for k in k_range
            data = full_data[k]
            write(io, string(k) *"\t"*string(data["lCI"])*"\t"*string(data["uCI"])*"\t"*string(data["mean"])*"\t"*string(data["median"])*"\n")
        end
    end
end

"""
    unresolveable_topological_incompatibilities(t1::Tree{T}, t2::Tree{T}, MCCs::Vector{Vector{String}})

note that a function can be compatible but not fully compatible - for example if it is only compatible
up to resolution
"""
function unresolveable_topological_incompatibilities(t1::Tree{T}, t2::Tree{T}, MCCs::Vector{Vector{String}}) where T 
    incompatiblilities = 0.0
    for mcc in MCCs
        is_comp =MTKTools.is_topo_compatible(t1, t2, mcc)
        if !is_comp
            incompatiblilities += length(mcc)
        end
    end
    return incompatiblilities/length(t1.lleaves)
end

"""
    topological_incompatibilities(t1::Tree{T}, t2::Tree{T}, MCCs::Vector{Vector{String}})

note that a function can be compatible but not fully compatible - for example if it is only compatible
up to resolution
"""
function topological_incompatibilities(t1::Tree{T}, t2::Tree{T}, MCCs::Vector{Vector{String}}) where T 
    incompatiblilities = 0.0
    for mcc in MCCs
        is_comp =MTKTools.is_full_topo_compatible(t1, t2, mcc)
        if !is_comp
            incompatiblilities += length(mcc)
        end
    end
    return incompatiblilities/length(t1.lleaves)
end

function run_one_sim(no_lineages::Int, rec_rate::Float64, strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    results = Dict()
    
    c = get_c(res, rec_rate; n=no_lineages, simtype)
    r = 10^rec_rate
    true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)

    for no_trees in k_range
        rand_order = sample(1:8, no_trees, replace = false)
        unresolved_trees = [copy(t) for t in true_trees[rand_order]]
        unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

        i_trees = [copy(t) for t in unresolved_trees]
        i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(;nMCMC=250, parallel=false, strict, rounds, final_no_resolve, pre_resolve))
        if rounds==1 && final_no_resolve && !pre_resolve
            @assert all([SplitList(i_t) == SplitList(t) for (i_t, t) in zip(i_trees, unresolved_trees)])
        end
        loc = sample(1:no_trees, 2, replace = false)
        names = [t.label for t in i_trees[loc]]

        percentage_incompatibilities = topological_incompatibilities(i_trees[loc][1], i_trees[loc][2], TreeKnit.get(i_MCCs, names...))
        percentage_unresolvable_incompatibilities = unresolveable_topological_incompatibilities(i_trees[loc][1], i_trees[loc][2], TreeKnit.get(i_MCCs, names...))

        if rounds==1 && final_no_resolve && !pre_resolve 
            @assert percentage_unresolvable_incompatibilities == 0.0
        end

        loc_consist = sample(1:no_trees, 3, replace = false)
        names = [t.label for t in i_trees[loc_consist]]
        MCCs = [TreeKnit.get(i_MCCs, [names[1], names[2]]...), TreeKnit.get(i_MCCs, [names[1], names[3]]...), TreeKnit.get(i_MCCs, [names[2], names[3]]...)]
        trees_ = [copy(t) for t in i_trees[loc_consist]]
        consistency_rate = MTKTools.consistency_rate(MCCs..., trees_)
        results[no_trees] = (percentage_incompatibilities, percentage_unresolvable_incompatibilities, sum(consistency_rate)/3)
    end

    return results
end

function run_MCC_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; strict=true, simtype=:flu, res=0.3, k_range=2:8, rounds=2, final_no_resolve=false, pre_resolve=false)
    percentage_incompatibilities = Dict{Int, Vector{Float32}}()
    percentage_unresolvable_incompatibilities = Dict{Int, Vector{Float32}}()
    consistency_rates = Dict{Int, Vector{Float32}}()

    percentage_incompatibilities_dict = Dict()
    percentage_unresolvable_incompatibilities_dict = Dict()
    consistency_rates_dict = Dict()
    
    sim_results = Dict()
    for i in 1:no_sim
        sim_results[i] = Dagger.@spawn run_one_sim(no_lineages, rec_rate, strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    end
    for no_trees in k_range
        percentage_incompatibilities[no_trees] = Float32[]
        percentage_unresolvable_incompatibilities[no_trees] = Float32[]
        consistency_rates[no_trees] = Float32[]
    end
    for i in 1:no_sim
        sim_result = fetch(sim_results[i])
        for no_trees in k_range
            push!(percentage_incompatibilities[no_trees], sim_result[no_trees][1])
            push!(percentage_unresolvable_incompatibilities[no_trees], sim_result[no_trees][2])
            push!(consistency_rates[no_trees], sim_result[no_trees][3])
        end
    end
    for no_trees in k_range
        for (p,v_dict) in zip([percentage_incompatibilities, percentage_unresolvable_incompatibilities, consistency_rates], [percentage_incompatibilities_dict, percentage_unresolvable_incompatibilities_dict, consistency_rates_dict])
            values = sort(p[no_trees])
            CI = [values[Int(round(no_sim*0.05))+1], values[Int(round(no_sim*0.95))]]
            dict = Dict("lCI" =>CI[1], "mean" => sum(values)/no_sim, "uCI" =>CI[2], "median" => values[Int(round(no_sim*0.5))] )

            v_dict[no_trees] = dict
        end
    end  
    return  (percentage_incompatibilities_dict, percentage_unresolvable_incompatibilities_dict, consistency_rates_dict)
end

function main()
    # n = 50
    # r = 0.0
    # o = "results_test"
    # num_sim = 10
    # simt = "flu"
    # res = 0.3
    # strict = true
    # rounds = 2
    # krange_ = true
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    num_sim = parsed_args["sim"]
    strict = parsed_args["strict"]
    simt = parsed_args["simtype"]
    res = parsed_args["res"]
    rounds = parsed_args["rounds"]
    krange_ = parsed_args["krange"]
    final_no_resolve = parsed_args["final-no-resolve"]
    pre_resolve = parsed_args["pre-resolve"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    if krange_
        k_range = 2:8
    else
        k_range = [4,8]
    end
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r simtype $simt")
    incomp, unres_incomp, consist = run_MCC_accuracy_simulations(num_sim, n, r; strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    write_output(o, incomp, "topo_incomp", r; k_range)
    write_output(o, unres_incomp, "topo_unres_incomp", r; k_range)
    write_output(o, consist, "consistency", r; k_range)
end

main()


