using Plots
using CSV
using DataFrames
using JSON3

json_string = read("../config.json", String)

params = JSON3.read(json_string)

pre_recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
#recomb_rate = [0]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
pre_x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

function add_to_dict_vec(dict_, file, key)
    if !haskey(dict_, key)
        dict_[key] = file[(file.k .== key), :].mean
    else
        append!(dict_[key], file[(file.k .== key), :].mean)
    end
end

for x_axis in ["scaled", "true"]
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
                    add_to_dict_vec(ic_mean_, ic_file, 0)
                    add_to_dict_vec(c_mean_, c_file, 0)
                    add_to_dict_vec(c_mean_, c_file, 1)
                    for k in 2:8
                        add_to_dict_vec(ic_mean_, ic_file, k)
                        add_to_dict_vec(c_mean_, c_file, k)
                        add_to_dict_vec(rf_mean_, rf_file, k)
                    end
                end
                if x_axis =="scaled"
                    x_axis_title = "scaled reassortment rate (ρ)"
                    recomb_rate = pre_recomb_rate
                    x_ticks = pre_x_ticks
                    xpad = 70  # adjust function of font size
                else
                    if sim== "kingman"
                        alpha = 1.0
                    else
                        alpha = 0.2
                    end
                    recomb_rate = ((1/10000)*(params["n"]/2)^alpha) .* pre_recomb_rate
                    x_ticks = ((1/10000)*(params["n"]/2)^alpha) .* pre_x_ticks
                    x_axis_title = "true reassortment rate (r)"
                    xpad = ((1/10000)*(params["n"]/2)^alpha)*70  # adjust function of font size   
                end
                p1 = plot(recomb_rate, c_mean_[0] .* 100, widen = false, label="rMCCs, %correct", ylabel="% correct new splits", xlabel=x_axis_title , linecolor="black", link=:xaxis, title="% Resolved Polytomies for Average Tree in k-Tree Sample", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6)
                plot!(recomb_rate, ic_mean_[0] .* 100, label="rMCCs, %incorrect", ylabel="% correct new splits", xlabel=x_axis_title , linecolor="black",linestyle=:dash)
                for no_trees in 2:8
                    plot!(recomb_rate, c_mean_[no_trees] .* 100, label="k="*string(no_trees)*", %correct", ylabel="% correct new splits", xlabel=x_axis_title , linecolor=l_color[no_trees-1])
                    plot!(recomb_rate, ic_mean_[no_trees] .* 100, label="k="*string(no_trees)*", %incorrect", ylabel="% correct new splits", xlabel=x_axis_title , linecolor=l_color[no_trees-1], linestyle=:dash)
                end
                vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                p2 = plot(recomb_rate, c_mean_[1] ./ params["n"], widen = false, label="average size MCCs", ylabel="average size MCCs", xlabel=x_axis_title , linecolor="red", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6, yaxis = :log10)
                vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                l = @layout [a{0.7h} ; b]
                p = plot(p1, p2, layout = l, link = :recomb_rate)
                savefig(p, "Plots/PercentageCorrectResolution_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
                p = plot(recomb_rate, rf_mean_[2], label="k=2, improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title , linecolor=l_color[1], title="Improvement in RF distance - Average Tree in k-Tree Sample", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                for no_trees in 3:8
                    plot!(recomb_rate, rf_mean_[no_trees], label="k="*string(no_trees)*", improvement RF distance", ylabel="RF distance to true tree tree, unresolved - infered", xlabel=x_axis_title, linecolor=l_color[no_trees-1])
                end
                savefig(p, "Plots/PercentageRF_improvement_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
            end
        end
    end
end
