using Plots
using CSV
using DataFrames
using JSON3
#Consistency/TopologicalCompatibility/
json_string = read("consistent_config.json", String)

params = JSON3.read(json_string)

recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkgreen", "lightgreen", "darkblue", "lightblue", "purple", "magenta", "red", "orange"]
x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for consistent in params["CONSISTENT"]
    for sim in params["SIMTYPE"]
        for res in params["RES"]
            for strict in params["STRICT"]
                tc_mean_ = Dict()
                tc_ur_mean_ = Dict()
                const_mean_ = Dict()
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
                                filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*consistent*"_"*fr*"/results_topo_incomp_"*string(rec)*".txt"
                                tc_file = CSV.read(filename, DataFrame)
                                filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*consistent*"_"*fr*"/results_topo_unres_incomp_"*string(rec)*".txt"
                                tc_ur_file = CSV.read(filename, DataFrame)
                                filename = "results/results_"*sim*"_"*res*"_"*strict*"_"*r*"_"*pr*"_"*consistent*"_"*fr*"/results_consistency_"*string(rec)*".txt"
                                const_file = CSV.read(filename, DataFrame)
                                for k in [4,8]
                                    if !haskey(tc_mean_, (r, k, pr_val, fr_val))
                                        tc_mean_[(r,k, pr_val, fr_val)] = tc_file[(tc_file.k .== k), :].mean
                                    else
                                        append!(tc_mean_[(r,k, pr_val, fr_val)], tc_file[(tc_file.k .== k), :].mean)
                                    end
                                    if !haskey(tc_ur_mean_, (r,k, pr_val, fr_val))
                                        tc_ur_mean_[(r,k, pr_val, fr_val)] = tc_ur_file[(tc_ur_file.k .== k), :].mean
                                    else
                                        append!(tc_ur_mean_[(r,k, pr_val, fr_val)], tc_ur_file[(tc_ur_file.k .== k), :].mean)
                                    end
                                    if !haskey(const_mean_, (r,k, pr_val, fr_val))
                                        const_mean_[(r,k, pr_val, fr_val)] = const_file[(const_file.k .== k), :].mean
                                    else
                                        append!(const_mean_[(r,k, pr_val, fr_val)], const_file[(const_file.k .== k), :].mean)
                                    end
                                end
                            end
                        end
                    end
                end
                
                for k in [4,8]
                    p = plot(recomb_rate, tc_mean_[(params["ROUNDS"][1], k, 0, 0)], label="r="*string(params["ROUNDS"][1])*", %incomp.", ylabel="% leaves in topo. incomp. MCC", xlabel="recombination rate", linecolor=l_color[1], title="% Topological Incompatibilities for Average Tree in "*string(k)*"-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)    
                    for r in params["ROUNDS"]
                        pos = parse(Int,r) + 4*(parse(Int,r)-1)
                        if r != params["ROUNDS"][1]
                            plot!(recomb_rate, tc_mean_[(r, k, 0, 0)], label="r="*string(r)*", %incomp.", ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                        end
                        plot!(recomb_rate, tc_ur_mean_[(r, k, 0, 0)], label="r="*string(r)*", %unres. incomp.", ylabel="", xlabel="recombination rate", linecolor=l_color[pos],linestyle=:dash)
                        for (i, (vals, name)) in enumerate(zip([(1,0), (0,1), (1,1)], ["pre-r", "final-no-r", "pre-r, final-no-r"]))
                            pos = parse(Int,r) +i+ 4*(parse(Int,r) -1)
                            plot!(recomb_rate, tc_mean_[(r, k, vals...)], label="r="*string(r)*", %incomp., "*name, ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                            plot!(recomb_rate, tc_ur_mean_[(r, k, vals...)], label="r="*string(r)*", %unres. incomp. "*name, ylabel="", xlabel="recombination rate", linecolor=l_color[pos],linestyle=:dash)
                        end
                        savefig(p, "Plots/Topo_Compat_k"*string(k)*"_"*sim*"_"*res*"_"*strict*"_"*consistent*".png")
                    end
                    p = plot(recomb_rate, const_mean_[(params["ROUNDS"][1], k, 0, 0)], label="r="*string(params["ROUNDS"][1])*", %incons. splits.", ylabel="%inconsistent splits", xlabel="recombination rate", linecolor=l_color[1], title="%MCC Inconsistencies for Average Tree in "*string(k)*"-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)    
                    for r in params["ROUNDS"]
                        pos = parse(Int,r) + 4*(parse(Int,r)-1)
                        if r != params["ROUNDS"][1]
                            plot!(recomb_rate, const_mean_[(r, k, 0, 0)], label="r="*string(r)*", %incons. splits", ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                        end
                        for (i, (vals, name)) in enumerate(zip([(1,0), (0,1), (1,1)], ["pre-r", "final-no-r", "pre-r, final-no-r"]))
                            pos = parse(Int,r) +i+ 4*(parse(Int,r) -1)
                            plot!(recomb_rate, const_mean_[(r, k, vals...)], label="r="*string(r)*", %incons. splits, "*name, ylabel="", xlabel="recombination rate", linecolor=l_color[pos])
                        end
                        savefig(p, "Plots/Consistency_k"*string(k)*"_"*sim*"_"*res*"_"*strict*"_"*consistent*".png")
                    end
                end
            end
        end
    end
end