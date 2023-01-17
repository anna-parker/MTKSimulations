using Plots
using CSV
using DataFrames
using JSON3

json_string = read("rounds_config.json", String)

params = JSON3.read(json_string)

pre_recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "darkgreen", "green", "lightgreen", "purple", "magenta", "red"]
pre_x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for x_axis in ["scaled", "true"]
    for sim in params["SIMTYPE"]
        for res in params["RES"]
            for strict in params["STRICT"]
                ic_mean_ = Dict()
                c_mean_ = Dict()
                rf_mean_ = Dict()
                for r in params["ROUNDS"]
                    for rec in params["REC"]
                        filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*string(r)*"/results_correct_"*string(rec)*".txt"
                        c_file = CSV.read(filename, DataFrame)
                        filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*string(r)*"/results_incorrect_"*string(rec)*".txt"
                        ic_file = CSV.read(filename, DataFrame)
                        filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*string(r)*"/results_rf_"*string(rec)*".txt"
                        rf_file = CSV.read(filename, DataFrame)
                        for k in 2:8
                            if !haskey(ic_mean_, (r,k))
                                ic_mean_[(r,k)] = ic_file[(ic_file.k .== k), :].mean
                            else
                                append!(ic_mean_[(r,k)], ic_file[(ic_file.k .== k), :].mean)
                            end
                            if !haskey(c_mean_, (r,k))
                                c_mean_[(r,k)] = c_file[(c_file.k .== k), :].mean
                            else
                                append!(c_mean_[(r,k)], c_file[(c_file.k .== k), :].mean)
                            end
                            if !haskey(rf_mean_, (r,k))
                                rf_mean_[(r,k)] = rf_file[(rf_file.k .== k), :].mean
                            else
                                append!(rf_mean_[(r,k)], rf_file[(rf_file.k .== k), :].mean)
                            end
                        end
                    end
                end

                if x_axis =="scaled"
                    x_axis_title = "scaled reassortment rate (ρ)"
                    recomb_rate = pre_recomb_rate
                    x_ticks = pre_x_ticks
                else
                    if sim== "kingman"
                        alpha = 1.0
                    else
                        alpha = 0.2
                    end
                    recomb_rate = ((1/10000)*(params["n"]/2)^alpha) .* pre_recomb_rate
                    x_ticks = ((1/10000)*(params["n"]/2)^alpha) .* pre_x_ticks
                    x_axis_title = "true reassortment rate (r)" 
                end
            
                p = plot(recomb_rate, c_mean_[(1,2)], label="k=2, 1 round, %correct", ylabel="% correct new splits", xlabel=x_axis_title, linecolor=l_color[1], title="% Resolved Polytomies for Average Tree in k-Tree Sample", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                plot!(recomb_rate, ic_mean_[(1,2)], label="k=2, 1 round, %incorrect", ylabel="", xlabel=x_axis_title, linecolor=l_color[1],linestyle=:dash)    
                for no_trees in [4,8]
                    pos = (Int(log2(no_trees))-1)*3 + 1
                    plot!(recomb_rate, c_mean_[(1, no_trees)], label="k="*string(no_trees)*", 1 round, %correct", ylabel="% correct new splits", xlabel=x_axis_title, linecolor=l_color[pos])
                    plot!(recomb_rate, ic_mean_[(1, no_trees)], label="k="*string(no_trees)*", 1 rounds, %incorrect", ylabel="", xlabel=x_axis_title, linecolor=l_color[pos], linestyle=:dash)
                end
                for r in params["ROUNDS"][2:end]
                    plot!(recomb_rate, c_mean_[(r, 2)], label="k=2, "*string(r)*" rounds, %correct", ylabel="% correct new splits", xlabel=x_axis_title, linecolor=l_color[r])
                    plot!(recomb_rate, ic_mean_[(r, 2)], label="k=2, "*string(r)*" rounds, %incorrect", ylabel="", xlabel=x_axis_title, linecolor=l_color[r], linestyle=:dash)
                    for no_trees in [4,8]
                        pos = (Int(log2(no_trees))-1)*3 + r
                        plot!(recomb_rate, c_mean_[(r, no_trees)], label="k="*string(no_trees)*", "*string(r)*" rounds, %correct", ylabel="% correct new splits", xlabel=x_axis_title, linecolor=l_color[pos])
                        plot!(recomb_rate, ic_mean_[(r, no_trees)], label="k="*string(no_trees)*", "*string(r)*" rounds, %incorrect", ylabel="", xlabel=x_axis_title, linecolor=l_color[pos], linestyle=:dash)
                    end
                end
                savefig(p, "Plots/PercentageCorrectResolution_"*sim*"_"*res*"_"*strict*".png")
                p = plot(recomb_rate, rf_mean_[(1, 2)], label="k=2, 1 round, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title, linecolor=l_color[1], title="Improvement in RF distance - Average Tree in k-Tree Sample", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                for no_trees in [4,8]
                    pos = (Int(log2(no_trees))-1)*3 + 1
                    plot!(recomb_rate, rf_mean_[(1, no_trees)], label="k="*string(no_trees)*", 1 round, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title, linecolor=l_color[pos])
                end
                for r in params["ROUNDS"][2:end]
                    plot!(recomb_rate, rf_mean_[(r, 2)], label="k=2, "*string(r)*" rounds, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title, linecolor=l_color[r])
                    for no_trees in [4,8]
                        pos = (Int(log2(no_trees))-1)*3 + r
                        plot!(recomb_rate, rf_mean_[(r, no_trees)], label="k="*string(no_trees)*", "*string(r)*" rounds, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title, linecolor=l_color[pos])
                    end
                end
                savefig(p, "Plots/PercentageRF_improvement_"*sim*"_"*res*"_"*strict*".png")
            end
        end
    end
end