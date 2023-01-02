using Plots
using CSV
using DataFrames
using JSON3

json_string = read("consistent_config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkgreen", "lightgreen", "darkblue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for consistent in params["CONSISTENT"]
    for sim in params["SIMTYPE"]
        for res in params["RES"]
            for strict in params["STRICT"]
                c_mean_ = Dict()
                ic_mean_ = Dict()
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
                                filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*consistent*"_"*fr*"/results_correct_"*string(rec)*".txt"
                                c_file = CSV.read(filename, DataFrame)
                                filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*consistent*"_"*fr*"/results_incorrect_"*string(rec)*".txt"
                                ic_file = CSV.read(filename, DataFrame)
                                for k in [2, 4,8]
                                    if !haskey(c_mean_, (r, k, pr_val, fr_val))
                                        c_mean_[(r,k, pr_val, fr_val)] = c_file[(c_file.k .== k), :].mean
                                    else
                                        append!(c_mean_[(r,k, pr_val, fr_val)], c_file[(c_file.k .== k), :].mean)
                                    end
                                    if !haskey(ic_mean_, (r,k, pr_val, fr_val))
                                        ic_mean_[(r,k, pr_val, fr_val)] = ic_file[(ic_file.k .== k), :].mean
                                    else
                                        append!(ic_mean_[(r,k, pr_val, fr_val)], ic_file[(ic_file.k .== k), :].mean)
                                    end
                                end
                            end
                        end
                    end
                end
                
                for k in [2,4,8]
                    p = plot(recomb_rate, c_mean_[(params["ROUNDS"][1], k, 0, 0)], label="r="*string(params["ROUNDS"][1])*", %correct", ylabel="%correct new splits", xlabel="recombination rate", linecolor=l_color[1], title="% Resolved Polytomies for Average Tree in "*string(k)*"-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)    
                    for r in params["ROUNDS"]
                        pos = 1 + 4*(parse(Int,r)-1)
                        if r != params["ROUNDS"][1]
                            plot!(recomb_rate, c_mean_[(r, k, 0, 0)], label="r="*string(r)*", %correct", ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                        end
                        plot!(recomb_rate, ic_mean_[(r, k, 0, 0)], label="r="*string(r)*", %incorrect", ylabel="", xlabel="recombination rate", linecolor=l_color[pos],linestyle=:dash)
                        for (i, (vals, name)) in enumerate(zip([(1,0), (0,1), (1,1)], ["pre-r", "final-no-r", "pre-r, final-no-r"]))
                            pos = 1 +i+ 4*(parse(Int,r) -1)
                            plot!(recomb_rate, c_mean_[(r, k, vals...)], label="r="*string(r)*", %correct, "*name, ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                            plot!(recomb_rate, ic_mean_[(r, k, vals...)], label="r="*string(r)*", %incorrect. "*name, ylabel="", xlabel="recombination rate", linecolor=l_color[pos],linestyle=:dash)
                        end
                        savefig(p, "Plots/PercentageCorrectResolution_k"*string(k)*"_"*sim*"_"*res*"_"*strict*"_"*consistent*".png")
                    end
                end
            end
        end
    end
end