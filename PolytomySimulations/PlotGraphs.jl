using Plots
using CSV
using DataFrames
using JSON3

json_string = read("../config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
#recomb_rate = [0]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for sim in params["SIMTYPE"]
    for res in params["RES"]
        for strict in params["STRICT"]
            ic_mean_ = Dict()
            c_mean_ = Dict()
            rf_mean_ = Dict()
            for rec in params["REC"]
                filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_correct_"*string(rec)*".txt"
                c_file = CSV.read(filename, DataFrame)
                filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_incorrect_"*string(rec)*".txt"
                ic_file = CSV.read(filename, DataFrame)
                filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_rf_"*string(rec)*".txt"
                rf_file = CSV.read(filename, DataFrame)
                for k in 2:8
                    if !haskey(ic_mean_, k)
                        ic_mean_[k] = ic_file[(ic_file.k .== k), :].mean
                    else
                        append!(ic_mean_[k], ic_file[(ic_file.k .== k), :].mean)
                    end
                    if !haskey(c_mean_, k)
                        c_mean_[k] = c_file[(c_file.k .== k), :].mean
                    else
                        append!(c_mean_[k], c_file[(c_file.k .== k), :].mean)
                    end
                    if !haskey(rf_mean_, k)
                        rf_mean_[k] = rf_file[(rf_file.k .== k), :].mean
                    else
                        append!(rf_mean_[k], rf_file[(rf_file.k .== k), :].mean)
                    end
                end
            end
        
            p = plot(recomb_rate, c_mean_[2], label="k=2, %correct", ylabel="% correct new splits", xlabel="recombination rate", linecolor=l_color[1], title="% Resolved Polytomies for Average Tree in k-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
            plot!(recomb_rate, ic_mean_[2], label="k=2, %incorrect", ylabel="", xlabel="recombination rate", linecolor=l_color[1],linestyle=:dash)
            for no_trees in 3:8
                plot!(recomb_rate, c_mean_[no_trees], label="k="*string(no_trees)*", %correct", ylabel="% correct new splits", xlabel="recombination rate", linecolor=l_color[no_trees-1])
                plot!(recomb_rate, ic_mean_[no_trees], label="k="*string(no_trees)*", %incorrect", ylabel="", xlabel="recombination rate", linecolor=l_color[no_trees-1], linestyle=:dash)
            end
            savefig(p, "Plots/PercentageCorrectResolution_"*sim*"_"*res*"_"*strict*".png")
            p = plot(recomb_rate, rf_mean_[2], label="k=2, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel="recombination rate", linecolor=l_color[1], title="Improvement in RF distance - Average Tree in k-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
            for no_trees in 3:8
                plot!(recomb_rate, rf_mean_[no_trees], label="k="*string(no_trees)*", improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel="recombination rate", linecolor=l_color[no_trees-1])
            end
            savefig(p, "Plots/PercentageRF_improvement_"*sim*"_"*res*"_"*strict*".png")
        end
    end
end