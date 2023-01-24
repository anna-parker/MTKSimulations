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

function write_output(outfolder, data_standard, data_size, type, rec_rate; k_range=2:8)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    filename = "results_"*type*"_"*string(rec_rate)*".txt"
    write_output_to_file(outfolder*"/"*filename, data_standard, data_size, ; k_range)
end

function write_output_to_file(filename, full_data, data_size; k_range=2:8)
    open(filename, "w") do io
        write(io, "k\tlCI\tuCI\tmean\tmedian\n")
        for k in k_range
            data = full_data[k]
            write(io, string(k) *"\t"*string(data["lCI"])*"\t"*string(data["uCI"])*"\t"*string(data["mean"])*"\t"*string(data["median"])*"\n")
        end
        new_range = vcat([1], collect(k_range))
        for k in new_range
            data = data_size[k]
            write(io, string(k+8) *"\t"*string(data["lCI"])*"\t"*string(data["uCI"])*"\t"*string(data["mean"])*"\t"*string(data["median"])*"\n")
        end

    end
end

function homogenity(a, b)
    leaves = sort(vcat(a...))
	assignments = [TestRecombTools.assignment_vector(leaves, mccs) for mccs in [a,b]]
    A = Clustering.counts(assignments[1], assignments[2])
    N = sum(A)
    (N == 0.0) && return 0.0

    entA = entropy(A)
    entArows = entropy(sum(A, dims=2))
    entAcols = entropy(sum(A, dims=1))

    hck = (entA - entAcols)/N
    hc = entArows/N + log(N)

    # Homogeneity
    h = hc == 0.0 ? 1.0 : 1.0 - hck/hc

    return h
end

function complete(a, b)
    leaves = sort(vcat(a...))
	assignments = [TestRecombTools.assignment_vector(leaves, mccs) for mccs in [a,b]]
    A = Clustering.counts(assignments[1], assignments[2])
    N = sum(A)
    (N == 0.0) && return 0.0

    entA = entropy(A)
    entArows = entropy(sum(A, dims=2))
    entAcols = entropy(sum(A, dims=1))

    hkc = (entA - entArows)/N
    hk = entAcols/N + log(N)

    # Completeness
    c = hk == 0.0 ? 1.0 : 1.0 - hkc/hk
    return c
end

function run_one_sim(no_lineages::Int, rec_rate::Float64, type, strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    a_index_i_vector = Dict()
    size_vector = Dict()
    
    c = get_c(res, rec_rate; n=no_lineages, simtype)
    r = 10^rec_rate
    true_trees, arg = MTKTools.get_trees(8, no_lineages; c, ρ = r, simtype)
    rMCCs = MTKTools.get_real_MCCs(8, arg)
    rMCCs = TreeKnit.MCC_set(8, [t.label for t in true_trees], rMCCs)
    for no_trees in k_range
        rand_order = sample(1:8, no_trees, replace = false)
        unresolved_trees = [copy(t) for t in true_trees[rand_order]]
        unresolved_trees = MTKTools.remove_branches(unresolved_trees; c)

        i_trees = [copy(t) for t in unresolved_trees]

        i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(;nMCMC=250, parallel=false, strict, rounds, final_no_resolve, pre_resolve))
        loc = sample(1:no_trees, 2, replace = false)
        names = [t.label for t in true_trees[rand_order][loc]]

        if type=="VI"
            a_index_i = TestRecombTools.varinfo_similarity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
        elseif type=="v-measure-complete" ##higher values means more similar and complete
            a_index_i = complete(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
        elseif type=="v-measure-homogenity"
            a_index_i = homogenity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
        else
            a_index_i = TestRecombTools.rand_index_similarity(TreeKnit.get(rMCCs, names...), TreeKnit.get(i_MCCs, names...))
        end

        a_index_i_vector[no_trees] = a_index_i
        MCCs = TreeKnit.get(i_MCCs, names...)
        size_infered_MCCs = sum([length(m) for m in MCCs])/length(MCCs)
        size_vector[no_trees] = size_infered_MCCs
    end
    rand_MCC = sample(1:8, 2, replace = false)
    MCCs = get(rMCCs, rand_MCC...)
    average_size_MCC = sum([length(m) for m in MCCs])/length(MCCs)
    size_vector[1] = average_size_MCC

    return a_index_i_vector, size_vector
end

function run_MCC_accuracy_simulations(no_sim::Int, no_lineages::Int, rec_rate::Float64; type="VI", strict=true, simtype=:flu, res=0.3, k_range=2:8, rounds=2, final_no_resolve=false, pre_resolve=false)
    average_accuracy_i = Dict()
    average_size_i = Dict()
    accuracy_index_i = Dict{Int, Vector{Float32}}()
    size_index_i = Dict{Int, Vector{Float32}}()
    sim_results = Dict()
    new_range = vcat([1], collect(k_range))
    for i in 1:no_sim
        sim_results[i] = Dagger.@spawn run_one_sim(no_lineages, rec_rate, type, strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    end
    size_index_i[1] = Float32[]
    for no_trees in k_range
        accuracy_index_i[no_trees] = Float32[]
        size_index_i[no_trees] = Float32[]
    end
    for i in 1:no_sim
        sim_result = fetch(sim_results[i])
        for no_trees in k_range
            push!(accuracy_index_i[no_trees], sim_result[1][no_trees])
        end
        for no_trees in new_range
            push!(size_index_i[no_trees], sim_result[2][no_trees])
        end
    end
    for no_trees in k_range
        accuracy_vector_i = sort(accuracy_index_i[no_trees])
        CI = [accuracy_vector_i[Int(round(no_sim*0.05))+1], accuracy_vector_i[Int(round(no_sim*0.95))]]
        dict = Dict("lCI" =>CI[1], "mean" => sum(accuracy_vector_i)/no_sim, "uCI" =>CI[2], "median" => accuracy_vector_i[Int(round(no_sim*0.5))] )

        average_accuracy_i[no_trees] = dict
    end  
    for no_trees in new_range
        size_vector_i = sort(size_index_i[no_trees])
        CI = [size_vector_i[Int(round(no_sim*0.05))+1], size_vector_i[Int(round(no_sim*0.95))]]
        dict = Dict("lCI" =>CI[1], "mean" => sum(size_vector_i)/no_sim, "uCI" =>CI[2], "median" => size_vector_i[Int(round(no_sim*0.5))] )

        average_size_i[no_trees] = dict
    end 
    return  average_accuracy_i, average_size_i
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
    # rounds = 2
    # krange_ = true
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
        k_range = [2,4,8]
    end
    println("Simulating ARGs and sequences of sample size $n and recombination rate $r simtype $simt")

    vi_accuracy_i, size_mccs_i = run_MCC_accuracy_simulations(num_sim, n, r; type=metric, strict, simtype, res, k_range, rounds, final_no_resolve, pre_resolve)
    write_output(o, vi_accuracy_i, size_mccs_i, metric, r; k_range)


end

main()


