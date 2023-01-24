using Plots
using CSV
using DataFrames
using JSON3

json_string = read("consistent_config.json", String)

params = JSON3.read(json_string)

pre_recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkgreen", "lightgreen", "darkblue", "lightblue", "purple", "magenta", "red", "orange"]
pre_x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for x_axis in ["scaled", "true"]
    for metric in params["METRIC"]
            for sim in params["SIMTYPE"]
                for res in params["RES"]
                    for strict in params["STRICT"]
                        tc_mean_ = Dict()
                        for rec in params["REC"]
                            for r in params["ROUNDS"]
                                for pr in params["PRE_RESOLVE"]
                                    for fr in params["FINAL_NO_RESOLVE"]
                                        if pr == "true"
                                            pr_val = 1
                                        else
                                            pr_val = 0
                                        end
                                        if fr == "true"
                                            fr_val = 1
                                        else
                                            fr_val = 0
                                        end
                                        filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*fr*"/results_"*metric*"_"*string(rec)*".txt"
                                        tc_file = CSV.read(filename, DataFrame)
                                        for k in [4,8]
                                            if !haskey(tc_mean_, (r, k, pr_val, fr_val))
                                                tc_mean_[(r,k, pr_val, fr_val)] = tc_file[(tc_file.k .== k), :].mean
                                            else
                                                append!(tc_mean_[(r,k, pr_val, fr_val)], tc_file[(tc_file.k .== k), :].mean)
                                            end
                                        end
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

                        if metric in ["v-measure-complete", "v-measure-homogenity"]
                            for k in [4,8]
                                p = plot(recomb_rate, tc_mean_[(params["ROUNDS"][1], k, 0, 0)], label="r="*string(params["ROUNDS"][1])*", "*metric, ylabel=metric, xlabel=x_axis_title, linecolor=l_color[1], title=metric*" infered to true MCCs in "*string(k)*"-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)    
                                for r in params["ROUNDS"]
                                    pos = 1 + 4*(parse(Int,r)-1)
                                    if r != params["ROUNDS"][1]
                                        plot!(recomb_rate, tc_mean_[(r, k, 0, 0)], label="r="*string(r)*", "*metric, ylabel="", xlabel=x_axis_title, linecolor=l_color[pos])
                                    end
                                    for (i, (vals, name)) in enumerate(zip([(1,0), (0,1), (1,1)], ["pre-r", "final-no-r", "pre-r, final-no-r"]))
                                        pos = 1 +i+ 4*(parse(Int,r) -1)
                                        plot!(recomb_rate, tc_mean_[(r, k, vals...)], label="r="*string(r)*", "*metric*", "*name, ylabel="", xlabel=x_axis_title, linecolor=l_color[pos])
                                    end
                                    savefig(p, "Plots/"*metric*"_k"*string(k)*"_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
                                end
                            end
                        else
                            for k in [4,8]
                                p = plot(recomb_rate, tc_mean_[(params["ROUNDS"][1], k, 0, 0)], label="r="*string(params["ROUNDS"][1])*", "*metric*" dist.", ylabel= metric*" distance", xlabel=x_axis_title, linecolor=l_color[1], title= metric*" dist infered from true MCCs in "*string(k)*"-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)    
                                for r in params["ROUNDS"]
                                    pos = 1 + 4*(parse(Int,r)-1)
                                    if r != params["ROUNDS"][1]
                                        plot!(recomb_rate, tc_mean_[(r, k, 0, 0)], label="r="*string(r)*", "*metric*" dist.", ylabel="", xlabel=x_axis_title, linecolor=l_color[pos])
                                    end
                                    for (i, (vals, name)) in enumerate(zip([(1,0), (0,1), (1,1)], ["pre-r", "final-no-r", "pre-r, final-no-r"]))
                                        pos = 1 +i+ 4*(parse(Int,r) -1)
                                        plot!(recomb_rate, tc_mean_[(r, k, vals...)], label="r="*string(r)*", "*metric*" dist., "*name, ylabel="", xlabel=x_axis_title, linecolor=l_color[pos])
                                    end
                                    savefig(p, "Plots/"*metric*"_k"*string(k)*"_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
                                end
                            end
                        end
                    end
                end
            end
        end
end
