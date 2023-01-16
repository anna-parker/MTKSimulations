using Plots
using CSV
using DataFrames
using JSON3

METRIC = ["rand", "VI", "v-measure-complete", "v-measure-homogenity"]

json_string = read("mcc_config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for sim in params["SIMTYPE"]
    for res in params["RES"]
        for strict in params["STRICT"]
            for metric in METRIC
                _mean = Dict()
                for rec in params["REC"]
                    filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_"*metric*"_"*string(rec)*".txt"
                    file = CSV.read(filename, DataFrame)
                    for k in 2:8
                        if !haskey(_mean, k)
                            _mean[k] = file[(file.k .== k), :].mean
                        else
                            append!(_mean[k], file[(file.k .== k), :].mean)
                        end
                    end
                end
                p = plot(recomb_rate, _mean[2], label="k=2, "*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[1], title=metric*" dist infered from true MCC in k-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                for no_trees in 3:8
                    plot!(recomb_rate, _mean[no_trees], label="k="*string(no_trees)*", "*metric*" dist.", ylabel=metric*" distance", xlabel="recombination rate", linecolor=l_color[no_trees-1])
                end
                savefig(p, "Plots/"*metric*"_Accuracy_"*sim*"_"*res*"_"*strict*".png")
            end
        end
    end
end



