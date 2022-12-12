using Plots
using CSV
using DataFrames
using JSON3

json_string = read("div_times_config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for sim in params["SIMTYPE"]
    for res in params["RES"]
        for strict in params["STRICT"]
            mean_ = Dict()
            median_ = Dict()
            var_ = Dict()
            for (n, rec) in enumerate(params["REC"])
                for k in range(2,8)
                    filename = "DivergenceTimeEstimations/results/results_"*sim*"_"*res*"_"*strict*"/tt_summary_"*string(n-1)*"_"*string(k)*".csv"
                    df = CSV.read(filename, DataFrame)
                    nr = nrow(df)
                    if k==2
                        if !haskey(mean_, 1)
                            mean_[1] = [sum(df.mean_old)/nr]
                            median_[1] = [sum(df.median_old)/nr]
                            var_[1] = [sum(df.var_old)/nr]
                        else
                            append!(mean_[1], sum(df.mean_old)/nr)
                            append!(median_[1], sum(df.median_old)/nr)
                            append!(var_[1], sum(df.var_old)/nr)
                        end
                    end
                    if !haskey(mean_, k)
                        mean_[k] = [sum(df.mean_arg)/nr]
                    else
                        append!(mean_[k], sum(df.mean_arg)/nr)
                    end
                    if !haskey(median_, k)
                        median_[k] = [sum(df.median_arg)/nr]
                    else
                        append!(median_[k], sum(df.median_arg)/nr)
                    end
                    if !haskey(var_, k)
                        var_[k] = [sum(df.var_arg)/nr]
                    else
                        append!(var_[k], sum(df.var_arg)/nr)
                    end
                end
            end
            for (name, vector) in zip(["mean", "median", "var"], [mean_, median_, var_])
                p = plot(recomb_rate, vector[1], label="k=1 (standard)", ylabel=name, xlabel="recombination rate", linecolor="brown", title="Difference Divergence Time Estimates ("*name*")", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                for no_trees in 2:8
                    plot!(recomb_rate, vector[no_trees], label="k="*string(no_trees), ylabel=name, xlabel="recombination rate", linecolor=l_color[no_trees-1])
                end
                savefig(p, "Plots/"*name*"_divergence_times_"*sim*"_"*res*"_"*strict*".png")
            end
        end
    end
end