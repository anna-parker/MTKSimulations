using TreeKnit
using TreeTools
using StatsBase
using ArgParse
using MTKTools


function write_output(outfolder, tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict, rec_rate)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    filenames = ["results_tp_"*string(rec_rate)*".txt", "results_fp_"*string(rec_rate)*".txt", "results_fn_"*string(rec_rate)*".txt", "results_tn_"*string(rec_rate)*".txt", "results_true_shared_"*string(rec_rate)*".txt"]
    for (filename, data_full) in zip(filenames, [tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict])
        write_output_to_file(outfolder*"/"*filename, data_full)
    end
end

function write_output_to_file(filename, full_data)
    open(filename, "w") do io
        write(io, "k\tlCI\tuCI\tmean\tmedian\n")
        for k in 2:8
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
            default = 0.0
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
    end

    return parse_args(s)
end

function get_true_number_shared(tree, MCC)
    shared_branches_ = MTKTools.map_shared_branches(MCC, tree)
    true_shared = 0
    for n in nodes(tree)
        if !isroot(n)
            if shared_branches_[n.label]
                true_shared += 1
            end
        end
    end
    return true_shared 
end


function run_shared_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; strict=true, simtype=:flu, res=0.3)
    r = 10^rec_rate
    tp_, fp_, fn_, tn_, true_shared_ = Dict(), Dict(), Dict(), Dict(), Dict()
    tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict = Dict(), Dict(), Dict(), Dict(), Dict()
    
    for i in 1:no_sim
        c = get_c(res, rec_rate; n=no_lineages, simtype)
        true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)
        rMCCs = MTKTools.get_real_MCCs(8, arg)
        rMCCs = TreeKnit.MCC_set(8, [t.label for t in true_trees], rMCCs)
        for no_trees in 2:8
            rand_order = sample(1:8, no_trees, replace = false)
            unresolved_trees = [copy(t) for t in true_trees[rand_order]]
            unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

            i_trees = [copy(t) for t in unresolved_trees]
            i_MCCs = run_treeknit!(i_trees, TreeKnit.OptArgs(;nMCMC=250, parallel=true, strict))
            
            loc = sample(1:no_trees, 2, replace = false)
            labels = [t.label for t in true_trees[rand_order][loc]]
            tp, fp, fn, tn = MTKTools.accuracy_shared_branches(i_trees[loc[1]], true_trees[rand_order][loc[1]], TreeKnit.get(i_MCCs, labels...), TreeKnit.get(rMCCs, labels...))
	        true_shared = get_true_number_shared(true_trees[rand_order][loc[1]], TreeKnit.get(rMCCs, labels...))
            
            for (vector, value) in zip([tp_, fp_, fn_, tn_ , true_shared_],[tp, fp, fn, tn, true_shared] )
                if !haskey(vector, no_trees)
                    vector[no_trees] = [value]
                else
                    append!(vector[no_trees], [value])
                end
            end
        end
    end

    for no_trees in 2:8
        for (vector, dict) in zip([tp_, fp_, fn_, tn_, true_shared_ ],[tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict] )
            change_vector = sort(vector[no_trees])
            CI = [change_vector[Int(round(no_sim*0.05))+1], change_vector[Int(round(no_sim*0.95))]]
            dict_ = Dict("lCI" =>CI[1], "mean" => sum(change_vector)/no_sim, "uCI" =>CI[2], "median" => change_vector[Int(round(no_sim*0.5))] )

            dict[no_trees] = dict_
        end
    end  
    return  tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict
end

function main()
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    num_sim = parsed_args["sim"]
    strict = parsed_args["strict"]
    simt = parsed_args["simtype"]
    res = parsed_args["res"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    #n = 10
    #r = -0.75
    #o = "results_test"
    #num_sim = 500
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r")
    tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict = run_shared_accuracy_simulations(num_sim, n, r; strict, simtype, res)
    if !isdir(o)
        mkdir(o)
    end
    write_output(o, tp_dict, fp_dict, fn_dict, tn_dict, true_shared_dict, r)
end

main()
