using MTKTools
using TreeTools
using TreeKnit
using TreeKnit.MTK
using TestRecombTools
using ArgParse
using StatsBase
using Clustering


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
        "--metric"
            help = "which accuracy measure - rand, VI, v-measure-complete, v-measure-homogenity"
            arg_type = String
            default = "VI"
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
    end

    return parse_args(s)
end

function write_output(outfolder, data_standard, type, rec_rate; k_range=2:8)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    filename = "results_standard_"*type*"_"*string(rec_rate)*".txt"
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

function v_measure_similarity(MCCs...; scale=true, β = 1)
	leaves = sort(vcat(first(MCCs)...))
	assignments = try
		[TestRecombTools.assignment_vector(leaves, mccs) for mccs in MCCs]
	catch err
		println(MCCs[1])
		println()
		println(MCCs[2])
		println()
		println(leaves)
		error(err)
	end
	out = 0
	Z = 0
	for i in 1:length(MCCs), j in (i+1):length(MCCs)
		out += Clustering.vmeasure(assignments[i], assignments[j]; β)
		Z += 1
	end
	if scale
		out /= log(length(leaves))
	end
	return out / Z
end

function run_MCC_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; type="VI", strict=true, simtype=:flu, res=0.3, k_range=2:8, rounds=2)
    r = 10^rec_rate
    average_accuracy_i = Dict()
    
    accuracy_index_i = Dict{Int, Vector{Float32}}()

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
            i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(;nMCMC=250, consistent = false, parallel=true, strict, rounds))
            loc = sample(1:no_trees, 2, replace = false)
            names = [t.label for t in true_trees[rand_order][loc]]

            if type=="VI"
                a_index_i = TestRecombTools.varinfo_similarity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
            elseif type=="v-measure-complete" ##higher values means more similar and complete
                a_index_i = v_measure_similarity(TreeKnit.get(i_MCCs, names...), TreeKnit.get(rMCCs, names...); β = 0)
            elseif type=="v-measure-homogenity"
                a_index_i = v_measure_similarity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...); β = 0)
            else
                a_index_i = TestRecombTools.rand_index_similarity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
            end
            
            if !haskey(accuracy_index_i, no_trees)
                accuracy_index_i[no_trees] = [a_index_i]
            else
                append!(accuracy_index_i[no_trees], a_index_i)
            end
        end
    end
    for no_trees in k_range
        accuracy_vector_i = sort(accuracy_index_i[no_trees])
        CI = [accuracy_vector_i[Int(round(no_sim*0.05))+1], accuracy_vector_i[Int(round(no_sim*0.95))]]
        dict = Dict("lCI" =>CI[1], "mean" => sum(accuracy_vector_i)/no_sim, "uCI" =>CI[2], "median" => accuracy_vector_i[Int(round(no_sim*0.5))] )

        average_accuracy_i[no_trees] = dict
    end  
    return  average_accuracy_i
end

function main()
    # n = 50
    # r = 0.0
    # o = "results_test"
    # num_sim = 10
    # metric = "VI"
    # simt = "flu"
    # res = 0.3
    # strict = true
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    num_sim = parsed_args["sim"]
    metric = parsed_args["metric"]
    strict = parsed_args["strict"]
    simt = parsed_args["simtype"]
    res = parsed_args["res"]
    rounds = parsed_args["rounds"]
    krange_ = parsed_args["krange"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    if krange_
        k_range = 2:8
    else
        k_range = [2,4,8]
    end
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r simtype $simt")
    if metric=="VI"
        vi_accuracy_i = run_MCC_accuracy_simulations(num_sim, n, r; type=metric, strict, simtype, res, k_range, rounds)
        write_output(o, vi_accuracy_i, metric, r; k_range)
    else
        rf_accuracy_i =  run_MCC_accuracy_simulations(num_sim, n, r; type=metric, strict, simtype, res, k_range, rounds)
        write_output(o, rf_accuracy_i, metric, r; k_range)
    end
end

main()


