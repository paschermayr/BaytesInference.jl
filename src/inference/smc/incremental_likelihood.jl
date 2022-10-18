################################################################################
"""
$(SIGNATURES)
Compute cumulative incremental likelihood from SMC Diagnostics given range 'effective_iter'.

# Examples
```julia
```

"""
function compute_cumℓincrement(diagnostics::Vector{T}, effective_iter::R) where {T<:BaytesSMC.SMCDiagnostics, R<:Union{StepRange{Int64, Int64}, UnitRange{Int64}}}
    return cumsum([diagnostics[chain].ℓincrement for chain in effective_iter])
end

################################################################################
"""
$(SIGNATURES)
Plot cumulative incremental likelihood against each other.

# Examples
```julia
```

"""
function plot_cumℓincrement(
    #A Vector for different cumulative log likelihood increments, see 'compute_cumℓincrement'
    cumℓincrement::Vector{T},
    #X axis Vector (with dates) to specify timeframe in plots
    dates = 1:length(cumℓincrement),
    #Modelnames
    modelnames::Vector{String} = [string("Model ", iter) for iter in eachindex(cumℓincrement)],
    benchmarkmodel::Int64 = 1;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {T}
    ArgCheck.@argcheck length(cumℓincrement) == length(modelnames)

    plot_score = plot(;
        layout=(1, 1),
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        size=plot_default_size,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
#=
#Plot cumulative log incremental likelihood over time for each model
    for iter in eachindex(cumℓincrement)
        plot!(dates, cumℓincrement[iter],
        color = param_color[iter], legend=:topleft,
        label=modelnames[iter],
        subplot=1,
        )
    end
=#
#Plot Log Bayes Factor of all Models against chosen benchmark
    ℓbayes = [cumℓincrement[iter] .- cumℓincrement[benchmarkmodel] for iter in eachindex(cumℓincrement)]
    for iter in eachindex(cumℓincrement)
        Plots.plot!(
            dates, ℓbayes[iter],
            label= modelnames[iter], legend=:topleft,
            ylabel= string("Cum. Log Bayes Factor \nDifference to ", modelnames[benchmarkmodel]),
            color = palette=Plots.palette(param_color, length(cumℓincrement)+1)[iter],
            subplot=1
        )
    end
    return plot_score
end

################################################################################
#export
export compute_cumℓincrement, plot_cumℓincrement
