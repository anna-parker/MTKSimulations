using Plots
using CSV
using DataFrames
using JSON3

json_string = read("rounds_config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "darkgreen", "green", "lightgreen", "purple", "magenta", "red"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]
metric = "VI"

for sim in params["SIMTYPE"]
    for res in params["RES"]
        for strict in params["STRICT"]
            _mean = Dict()
            for r in params["ROUNDS"]
                for rec in params["REC"]
                    filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*string(r)*"/results_standard_"*metric*"_"*string(rec)*".txt"
                    file = CSV.read(filename, DataFrame)
                    for k in 2:8
                        if !haskey(_mean, (r, k))
                            _mean[(r, k)] = file[(file.k .== k), :].mean
                        else
                            append!(_mean[(r, k)], file[(file.k .== k), :].mean)
                        end
                    end
                end
            end
            p = plot(recomb_rate, _mean[(1,2)], label="k=2, 1 round"*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[1], title=metric*" dist infered from true MCC in k-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
            for no_trees in [4,8]
                pos = (Int(log2(no_trees))-1)*3 + 1
                plot!(recomb_rate, _mean[(1, no_trees)], label="k="*string(no_trees)*", 1 round "*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[pos])
                end
            for r in params["ROUNDS"][2:end]
                plot!(recomb_rate, _mean[(r, 2)], label="k=2, "*string(r)*" rounds"*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[r])
                for no_trees in [4,8]
                    pos = (Int(log2(no_trees))-1)*3 + r
                    plot!(recomb_rate, _mean[(r, no_trees)], label="k="*string(no_trees)*", "*string(r)*" rounds"*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[pos])
                end
            end
            savefig(p, "Plots/"*metric*"_Accuracy_"*sim*"_"*res*"_"*strict*".png")
        end
    end
end