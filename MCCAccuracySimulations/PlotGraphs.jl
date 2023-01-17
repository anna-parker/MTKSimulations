using Plots
using CSV
using DataFrames
using JSON3

METRIC = ["rand", "VI", "v-measure-complete", "v-measure-homogenity"]

json_string = read("mcc_config.json", String)

params = JSON3.read(json_string)

pre_recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange"]
pre_x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

for x_axis in ["scaled", "true"]
    for sim in params["SIMTYPE"]
        for res in params["RES"]
            for strict in params["STRICT"]
                for metric in METRIC
                    _mean_ = Dict()
                    for rec in params["REC"]
                        filename = "results/results_"*sim*"_"*res*"_"*strict*"/results_"*metric*"_"*string(rec)*".txt"
                        file = CSV.read(filename, DataFrame)
                        for k in 2:16
                            if !haskey(_mean_, k)
                                _mean_[k] = file[(file.k .== k), :].mean
                            else
                                append!(_mean_[k], file[(file.k .== k), :].mean)
                            end
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
                    if metric =="v-measure-complete"
                        label_name = "complete"
                    else if metric =="v-measure-homogenity"
                        label_name = "homogenity"
                    else
                        label_name = metric*" dist."
                    end
                    p1 = plot(recomb_rate, _mean_[2], widen = false, label="k=2, "*label_name, ylabel=label_name, xlabel=x_axis_title , linecolor=l_color[1], link=:xaxis, title=metric*" dist infered from true MCC in k-Tree Sample", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6)
                    for no_trees in 3:8
                        plot!(recomb_rate, _mean_[no_trees], label="k="*string(no_trees)*", "*label_name, ylabel=label_name, xlabel=x_axis_title , linecolor=l_color[no_trees-1])
                    end
                    vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                    
                    p2 = plot(recomb_rate, _mean_[9] ./ params["n"], widen = false, label="size real MCCs", ylabel="average size MCCs", xlabel=x_axis_title , linecolor="black", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6, yaxis = :log10)
                    for no_trees in 10:16
                        plot!(recomb_rate, _mean_[no_trees], label="size MCCs", ylabel="average size MCCs", xlabel=x_axis_title, linecolor=l_color[no_trees-9])
                    end
                    vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                    l = @layout [a{0.7h} ; b]
                    p = plot(p1, p2, layout = l, link = :recomb_rate)

                    savefig(p, "Plots/"*metric*"_Accuracy_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
                end
            end
        end
    end
end



