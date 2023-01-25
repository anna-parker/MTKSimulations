using TreeKnit
using TreeTools
using ArgParse
using StatsBase 
using MTKTools

function write_output(outfolder, data_correct, data_incorrect, rec_rate; k_range=2:8)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    filename_correct = "results_correct_"*string(rec_rate)*".txt"
    filename_incorrect = "results_incorrect_"*string(rec_rate)*".txt"
    for (filename, data_full) in zip([filename_correct, filename_incorrect], [data_correct, data_incorrect])
        write_output_to_file(outfolder*"/"*filename, data_full; k_range)
    end
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
        "--metric"
            help = "use RF distance as metric"
            arg_type = String
            default = "splits"
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
            default = true
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

function new_split_accuracy(true_trees, unres_trees, infered_trees)
    @assert length(true_trees) == length(infered_trees) == length(unres_trees)
    no_trees = length(true_trees)
    correct_new_splits_trees = []
    incorrect_new_splits_trees = []
    for i in 1:no_trees
        true_splits = SplitList(true_trees[i])
        unres_splits = SplitList(unres_trees[i])
        test_splits = SplitList(infered_trees[i])

        no_true_splits = length(true_splits)
        no_unres_splits = length(unres_splits)
        missing_splits = no_true_splits - no_unres_splits + 1
        @assert missing_splits >0 
        correct_new_splits = 0
        incorrect_new_splits = 0
        for s in test_splits
            if s ∈ true_splits && s ∉ unres_splits
                correct_new_splits +=1
            elseif s ∉ true_splits && s ∉ unres_splits
                incorrect_new_splits +=1 
            end
        end
        append!(correct_new_splits_trees, correct_new_splits/missing_splits)
        append!(incorrect_new_splits_trees, incorrect_new_splits/missing_splits)
    end
    return correct_new_splits_trees, incorrect_new_splits_trees
end

function new_tree_RF_dist(true_trees, test_trees)
    @assert length(true_trees) == length(test_trees)
    no_trees = length(true_trees)
    rf = []
    for i in 1:no_trees
        append!(rf, TreeTools.RF_distance(true_trees[i], test_trees[i]))
    end
    return rf
end

function resolve_using_realMCCs!(trees, rMCCs; strict=true)
    l_t = length(trees)
    for k in 1:8
        for i in 1:(l_t-1), j in (i+1):l_t
            MCCs = get(rMCCs, trees[i].label, trees[j].label)
            TreeKnit.resolve!(trees[i], trees[j], MCCs; tau = 0., strict)
        end
    end
end

function add_to_dict_vec(dict_, key, value)
    if !haskey(dict_, key)
        dict_[key] = [value]
    else
        append!(dict_[key], value)
    end
end

function run_tree_poly_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; simtype = :flu, res=0.3, strict=true, k_range=2:8, rounds=2, final_no_resolve=false, pre_resolve=false)
    r = 10^rec_rate
    percent_correct_new_splits_vector = Dict()
    percent_incorrect_new_splits_vector = Dict()

    percent_correct_new_splits_i = Dict{Int, Vector{Float32}}()
    percent_incorrect_new_splits_i = Dict{Int, Vector{Float32}}()
    
    for i in 1:no_sim
        c = get_c(res, rec_rate; n=no_lineages, simtype)
        true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)
        rMCCs = MTKTools.get_real_MCCs(8, arg)
        rMCCs = TreeKnit.MCC_set(8, [t.label for t in true_trees], rMCCs)
        for no_trees in k_range
            rand_order = sample(1:8, no_trees, replace = false)
            unresolved_trees = [copy(t) for t in true_trees[rand_order]]
            unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

            i_trees = [copy(t) for t in unresolved_trees]
            i_MCCs = run_treeknit!(i_trees, TreeKnit.OptArgs(;nMCMC=250, parallel=true, strict, rounds, final_no_resolve, pre_resolve))
            
            loc = rand(1:no_trees)
            correct_new_splits_i, incorrect_new_splits_i = new_split_accuracy([true_trees[rand_order][loc]], [unresolved_trees[loc]], [i_trees[loc]])
            
            add_to_dict_vec(percent_correct_new_splits_i, no_trees, correct_new_splits_i[1])
            add_to_dict_vec(percent_incorrect_new_splits_i, no_trees, incorrect_new_splits_i[1])

            if no_trees == 8
                i_trees_for_rMCCs_strict = [copy(t) for t in unresolved_trees]
                resolve_using_realMCCs!(i_trees_for_rMCCs_strict, rMCCs; strict=true)
                correct_new_splits_rMCC_strict, incorrect_new_splits_rMCC_strict = new_split_accuracy([true_trees[rand_order][loc]], [unresolved_trees[loc]], [i_trees_for_rMCCs_strict[loc]])
                add_to_dict_vec(percent_correct_new_splits_i, 0, correct_new_splits_rMCC_strict[1])
                add_to_dict_vec(percent_incorrect_new_splits_i, 0, incorrect_new_splits_rMCC_strict[1])

                rand_MCC = sample(1:8, 2, replace = false)
                MCCs = get(rMCCs, rand_MCC...)
                average_size_MCC = sum([length(m) for m in MCCs])/length(MCCs)
                add_to_dict_vec(percent_correct_new_splits_i, 1, average_size_MCC)
                add_to_dict_vec(percent_incorrect_new_splits_i, 1, average_size_MCC)
            end
        end
    end
    CI_new_splits_correct = []
    CI_new_splits_incorrect = []
    for no_trees in vcat([0, 1], k_range)
        cs_ = sort(percent_correct_new_splits_i[no_trees])
        ics_ = sort(percent_incorrect_new_splits_i[no_trees])
        CI_new_splits_correct = [cs_[Int(round(no_sim*0.05))+1], cs_[Int(round(no_sim*0.95))]]
        CI_new_splits_incorrect = (ics_[Int(round(no_sim*0.05))+1], ics_[Int(round(no_sim*0.95))])
        correct_data = Dict("lCI" =>CI_new_splits_correct[1], "mean" => sum(percent_correct_new_splits_i[no_trees])/no_sim, "uCI" =>CI_new_splits_correct[2], "median" => cs_[Int(round(no_sim*0.5))] )
        incorrect_data = Dict("lCI" =>CI_new_splits_incorrect[1], "mean" => sum(percent_incorrect_new_splits_i[no_trees])/no_sim, "uCI" =>CI_new_splits_incorrect[2], "median" => ics_[Int(round(no_sim*0.5))] )
        
        percent_correct_new_splits_vector[no_trees] = correct_data
        percent_incorrect_new_splits_vector[no_trees] = incorrect_data
    end  
    return  percent_correct_new_splits_vector, percent_incorrect_new_splits_vector
end

function run_tree_rf_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; simtype = :flu, res=0.3, strict=true, k_range=2:8, rounds=2, final_no_resolve=false, pre_resolve=false)
    r = 10^rec_rate
    rf_change_info = Dict()
    rf_change_vector = Dict()
    
    for i in 1:no_sim
        c = get_c(res, rec_rate; n=no_lineages, simtype)
        true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)
        for no_trees in k_range
            rand_order = sample(1:8, no_trees, replace = false)
            unresolved_trees = [copy(t) for t in true_trees[rand_order]]
            unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

            i_trees = [copy(t) for t in unresolved_trees]
            i_MCCs = run_treeknit!(i_trees, TreeKnit.OptArgs(; parallel=true, strict, rounds, final_no_resolve, pre_resolve))
            
            loc = rand(1:no_trees)
            old_rf = new_tree_RF_dist([true_trees[rand_order][loc]], [unresolved_trees[loc]])
            new_rf = new_tree_RF_dist([true_trees[rand_order][loc]], [i_trees[loc]])

            if old_rf[1] != 0
                change = (old_rf[1] - new_rf[1])/old_rf[1]
            else
                change = (old_rf[1] - new_rf[1])
            end

            add_to_dict_vec(rf_change_vector, no_trees, change)
        end
    end

    for no_trees in k_range
        change_vector = sort(rf_change_vector[no_trees])
        CI = [change_vector[Int(round(no_sim*0.05))+1], change_vector[Int(round(no_sim*0.95))]]
        dict = Dict("lCI" =>CI[1], "mean" => sum(change_vector)/no_sim, "uCI" =>CI[2], "median" => change_vector[Int(round(no_sim*0.5))] )

        rf_change_info[no_trees] = dict
    end  
    return  rf_change_info
end

function main()
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    num_sim = parsed_args["sim"]
    metric = parsed_args["metric"]
    simt = parsed_args["simtype"]
    res = parsed_args["res"]
    strict = parsed_args["strict"]
    rounds = parsed_args["rounds"]
    k_range_ = parsed_args["krange"]
    final_no_resolve = parsed_args["final-no-resolve"]
    pre_resolve = parsed_args["pre-resolve"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    if k_range_
        k_range = 2:8
    else
        k_range = [2,4,8]
    end
    # n = 10
    # r = 0.0
    # o = "results_test"
    # num_sim = 10
    # metric = "splits"
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r and k_range in $k_range")
    if metric=="rf"
        rf_info = run_tree_rf_accuracy_simulations(num_sim, n, r; simtype, res, strict, rounds, k_range, final_no_resolve, pre_resolve)
        if !isdir(o)
            mkdir(o)
        end
        write_output_to_file(o*"/results_rf_"*string(r)*".txt", rf_info; k_range)
    else
        percent_correct_new_splits_vector, percent_incorrect_new_splits_vector = run_tree_poly_accuracy_simulations(num_sim, n, r; simtype, res, strict, rounds, k_range, final_no_resolve, pre_resolve)
        new_range = vcat([0, 1], collect(k_range))
        write_output(o, percent_correct_new_splits_vector, percent_incorrect_new_splits_vector, r; k_range=new_range)
    end
end

main()
