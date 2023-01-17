using Plots
using CSV
using DataFrames
using JSON3

json_string = read("div_times_config.json", String)

params = JSON3.read(json_string)

pre_recomb_rate = [round(10^e, digits=4) for e in range(-4, 0, length(params["REC"]))]
l_color = ["darkblue", "blue", "lightblue", "purple", "magenta", "red", "orange", "brown"]
pre_x_ticks= [round(10^e, digits=4) for e in range(-4, 0, 5)]

function append_vector_to_dict(dict_, key, vector; length_=false)
    mean_vector = filter(p->!ismissing(p), vector)
    if length(mean_vector) == 0
        sum_ = NaN
    else
        sum_ = sum(mean_vector)/length(mean_vector)
    end
    if length_
        sum_ = (sum_/(2*params["n"] - 4))*100
    end
    if !haskey(dict_, key)
        dict_[key] = [sum_] 
    else
        append!(dict_[key], sum_ )
    end
end

for x_axis in ["scaled", "true"]
    for sim in params["SIMTYPE"]
        for res in params["RES"]
            for strict in params["STRICT"]
                ## initialize dictionaries
                mean_ = Dict()
                median_ = Dict()
                var_ = Dict()
                length_ = Dict()
                if params["SEGMENTS"] ==2
                    for val in ["TP", "TN", "FN", "FP"]
                        mean_[val] = Dict()
                        median_[val] = Dict()
                        var_[val] = Dict()
                        length_[val] = Dict()
                    end
                end
                for (n, rec) in enumerate(params["REC"])
                    for k in range(2,params["SEGMENTS"])
                        filename = "results/results_"*sim*"_"*res*"_"*strict*"/tt_summary_"*string(n-1)*"_"*string(k)*".csv"
                        df = CSV.read(filename, DataFrame)
                        nr = nrow(df)
                        if params["SEGMENTS"] == 2
                            for val in ["TP", "TN", "FN", "FP"]
                                append_vector_to_dict(mean_[val], 1, df[:, val*"_mean_old"])
                                append_vector_to_dict(median_[val], 1, df[:, val*"_median_old"])
                                append_vector_to_dict(var_[val], 1, df[:, val*"_var_old"])
                                append_vector_to_dict(length_[val], 1, df[:, val*"_length"]; length_=true)
                                append_vector_to_dict(mean_[val], k, df[:, val*"_mean_arg"])
                                append_vector_to_dict(median_[val], k, df[:, val*"_median_arg"])
                                append_vector_to_dict(var_[val], k, df[:, val*"_var_arg"])
                                append_vector_to_dict(length_[val], k, df[:, val*"_length"]; length_=true)
                            end
                        else 
                            if k==2 
                                append_vector_to_dict(mean_, 1, df.mean_old)
                                append_vector_to_dict(median_, 1, df.median_old)
                                append_vector_to_dict(var_, 1, df.var_old)
                            end
                            append_vector_to_dict(mean_, k, df.mean_arg)
                            append_vector_to_dict(median_, k, df.median_arg)
                            append_vector_to_dict(var_, k, df.var_arg)
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
                if params["SEGMENTS"] == 2
                    for (name, vector) in zip(["mean", "median", "var"], [mean_, median_, var_])
                        xpad = 20   # adjust function of font size
                        p1 = plot(recomb_rate, vector["TP"][1], widen = false, label="k=1, TP", ylabel=name, xlabel=x_axis_title, linecolor=l_color[2], linestyle=:dash, link=:xaxis!, title="Difference Divergence Time Estimates ("*name*")", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6)
                        for (i, val) in enumerate(["TP", "TN", "FN", "FP"])
                            pos = 2*i
                            if val !="TP"
                                plot!(recomb_rate, vector[val][1], label="k=1, "*val, ylabel=name, xlabel=x_axis_title, linecolor=l_color[pos], linestyle=:dash)
                            end
                            for no_trees in range(2,params["SEGMENTS"])
                                plot!(recomb_rate, vector[val][no_trees], label="k="*string(no_trees)*", "*val, ylabel=name, xlabel=x_axis_title, linecolor=l_color[pos])
                            end
                        end
                        vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                        p2 = plot(recomb_rate, length_["TP"][1], widen = false,  ylabel = "% shared branches", ylim = [0,100], xaxis= :log10, label="%TP", linecolor=l_color[2], link=:xaxis, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :topright, legendfontsize=6)
                        for (i, val) in enumerate(["TN", "FN", "FP"])
                            pos = 2*(i+1) 
                            plot!(recomb_rate, length_[val][1], label="%"*val,  ylabel = "% shared branches", xlabel=x_axis_title, linecolor=l_color[pos], linestyle=:dash)
                        end
                        vspan!([recomb_rate[end], recomb_rate[end]+xpad], c=:white, lc=:white, label=false)
                        l = @layout [a{0.7h} ; b]
                        p = plot(p1, p2, layout = l, link = :recomb_rate)
                        savefig(p, "Plots/"*name*"_divergence_times_"*sim*"_"*res*"_"*strict*".png")
                    end
                else
                    for (name, vector) in zip(["mean", "median", "var"], [mean_, median_, var_])
                        p = plot(recomb_rate, vector[1], label="k=1 (standard)", ylabel=name, xlabel=x_axis_title, linecolor="brown", title="Difference Divergence Time Estimates ("*name*")", titlefontsize=10, margin=8Plots.mm, xaxis= :log10, xguidefontsize=7, yguidefontsize=7, xtickfontsize=6, ytickfontsize=6, xticks=x_ticks, legend = :outertopleft, legendfontsize=6)
                        for no_trees in range(2,params["SEGMENTS"])
                            plot!(recomb_rate, vector[no_trees], label="k="*string(no_trees), ylabel=name, xlabel=x_axis_title, linecolor=l_color[no_trees-1])
                        end
                        savefig(p, "Plots/"*name*"_divergence_times_"*sim*"_"*res*"_"*strict*"_"*x_axis*".png")
                    end
                end
            end
        end
    end
end