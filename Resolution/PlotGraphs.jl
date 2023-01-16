using Plots
using CSV
using DataFrames
using JSON3

json_string = read("../test_config.json", String)

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
            mean_ = Dict()
            for rec in params["REC"]
                filename = "results/results_"*sim*"_"*res*"/results_"*string(rec)*".txt"
                c_file = CSV.read(filename, DataFrame)
                for k in 0:2
                    add_to_dict_vec(mean_, c_file, k)
                end
            end
            if x_axis =="scaled"
                x_axis_title = "scaled recombination rate (ρ)"
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
                x_axis_title = "true recombination rate (r)"
                xpad = ((1/10000)*(params["n"]/2)^alpha)*70  # adjust function of font size   
            end
            p1 = plot(recomb_rate, mean_[1], widen = false, label="rMCCs - # ambig. splits", ylabel="# new splits", xlabel=x_axis_title , linecolor="red", link=:xaxis, title="Average #Splits from real MCCs in Average Tree in 8-Tree ARG", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6)
            plot!(recomb_rate, mean_[2], label="rMCCs - # unambig. splits", ylabel="# new splits", xlabel=x_axis_title , linecolor="blue")

            vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
            p2 = plot(recomb_rate, mean_[0] ./ params["n"], widen = false, label="average size MCCs", ylabel="average size MCCs", xlabel=x_axis_title , linecolor="black", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6, yaxis = :log10)
            vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
            l = @layout [a{0.7h} ; b]
            p = plot(p1, p2, layout = l, link = :recomb_rate)
            savefig(p, "Plots/Resolution_"*sim*"_"*res*"_"*x_axis*".png")
        end
    end
end
