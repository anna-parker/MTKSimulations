using Plots
using CSV
using DataFrames
using JSON3

json_string = read("../config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]
total_no_branches = 2*params["n"] -2

for sim in params["SIMTYPE"]
    for res in params["RES"]
        for strict in params["STRICT"]
            tn_mean, fn_mean, tp_mean, fp_mean, true_shared_mean = Dict(), Dict(), Dict(), Dict(), Dict()
            for rec in params["REC"]
                for (type, mean_) in zip(["tn", "fn", "tp", "fp", "true_shared"], [tn_mean, fn_mean, tp_mean, fp_mean, true_shared_mean])
                    filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_"*type*"_"*string(rec)*".txt"
                    file = CSV.read(filename, DataFrame)
                    for k in 2:8
                        if !haskey(mean_, k)
                            mean_[k] = file[(file.k .== k), :].mean
                        else
                            append!(mean_[k], file[(file.k .== k), :].mean)
                        end
                    end
                end
            end

            p = plot(recomb_rate, (total_no_branches .- true_shared_mean[2])./ total_no_branches, label="true % not shared branches", ylabel="not shared branches", xlabel="recombination rate", linecolor="black", linestyle=:dash, title="Percentage not shared branches out of total branches", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
            plot!(recomb_rate, tn_mean[2]./ total_no_branches, label="k=2, % TP inferred", ylabel="not shared branches", xlabel="recombination rate", linecolor="darkblue")
            plot!(recomb_rate, fn_mean[2]./ total_no_branches, label="k=2, % FP", ylabel="not shared branches", xlabel="recombination rate", linecolor="darkblue", linestyle=:dash)
            for no_trees in 3:8
                plot!(recomb_rate, tn_mean[no_trees]./ total_no_branches, label="k="*string(no_trees)*", % TP", ylabel="not shared branches", xlabel="recombination rate", linecolor=l_color[no_trees-1])
                plot!(recomb_rate, fn_mean[no_trees]./ total_no_branches, label="k="*string(no_trees)*", % FP", ylabel="not shared branches", xlabel="recombination rate", linecolor=l_color[no_trees-1], linestyle=:dash)
            end
            savefig(p, "Plots/Percentage_not_shared_branches_"*sim*"_"*res*"_"*strict*".png")

            p = plot(recomb_rate, true_shared_mean[2]./ total_no_branches, label="true % shared branches", ylabel="shared branches", xlabel="recombination rate", linecolor="black", linestyle=:dash, title="Percentage shared branches out of total branches", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
            plot!(recomb_rate, tp_mean[2]./ total_no_branches, label="k=2, % TP inferred", ylabel="shared branches", xlabel="recombination rate", linecolor="darkblue")
            plot!(recomb_rate, fp_mean[2]./ total_no_branches, label="k=2, % FP", ylabel="shared branches", xlabel="recombination rate", linecolor="darkblue", linestyle=:dash)
            for no_trees in 3:8
                plot!(recomb_rate, tp_mean[no_trees]./ total_no_branches, label="k="*string(no_trees)*", % TP", ylabel="shared branches", xlabel="recombination rate", linecolor=l_color[no_trees-1])
                plot!(recomb_rate, fp_mean[no_trees]./ total_no_branches, label="k="*string(no_trees)*", % FP", ylabel="shared branches", xlabel="recombination rate", linecolor=l_color[no_trees-1], linestyle=:dash)
            end
            savefig(p, "Plots/Percentage_shared_branches_"*sim*"_"*res*"_"*strict*".png")
        end
    end
end

