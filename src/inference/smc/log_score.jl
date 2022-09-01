################################################################################
"""
$(SIGNATURES)
Compute log score of trace given data.

# Examples
```julia
```

"""
function compute_ℓscore(
    predictionᵛ::Vector{Vector{S}}, ℓweightsₙ::Vector{Vector{T}}, data::D; burnin=0
) where {S,T<:Real,D}
    if length(predictionᵛ[1][1]) > 1
        return println("Log score only supported for univariate data, use escore instead")
    end
    ## Initiate container
    #!NOTE: +2 because we start with ( size(latent, 1) - length(predictionᵛ) ) for warmup, then predict for +2 at +1 after warmup
    start_date = size(data, 1) - length(predictionᵛ) + burnin + 2
    #!NOTE: can only go up to T-1, as predictions at endpoint are for T+1, which we do not have.
    ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ = zeros(Float64, (length(predictionᵛ) - burnin - 1))
    ℓscoreʷᵉⁱᵍʰᵗᵉᵈ = zeros(Float64, (length(predictionᵛ) - burnin - 1))
    ## Loop through predictions
    #!NOTE: can only go up to T-1, as predictions at endpoint are for T+1, which we do not have.
    for iter in Base.OneTo((length(predictionᵛ) - burnin - 1))
        ## 1 Forecast to t+1
        predictionₜ₊₁ = predictionᵛ[iter + burnin]
        weightₜ₊₁ = exp.(ℓweightsₙ[iter + burnin])
        ## 2 Obtain weighted KDE from observed data
        kernelʷᵉⁱᵍʰᵗᵉᵈ = KernelDensity.kde(predictionₜ₊₁; weights=weightₜ₊₁, kernel=Normal)
        kernelᵘⁿʷᵉⁱᵍʰᵗᵉᵈ = KernelDensity.kde(predictionₜ₊₁; kernel=Normal)
        ## 3 Calculate Log scores
        ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ[iter] =
            -1 * log(KernelDensity.pdf(kernelʷᵉⁱᵍʰᵗᵉᵈ, data[start_date + iter - 1]))
        ℓscoreʷᵉⁱᵍʰᵗᵉᵈ[iter] =
            -1 * log(KernelDensity.pdf(kernelᵘⁿʷᵉⁱᵍʰᵗᵉᵈ, data[start_date + iter - 1]))
    end
    ## Return log scores
    return ℓscoreʷᵉⁱᵍʰᵗᵉᵈ, ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ
end

################################################################################
"""
$(SIGNATURES)
Plot log score of trace given data.

# Examples
```julia
```

"""
function plot_logscore(
    trace::Trace,
    smc,
    data::D;
    burnin=0,
    delete∞=false,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {D}
    ## Get predictions
    predictionᵛ = [
        trace.diagnostics[chain].base.prediction for chain in eachindex(trace.diagnostics)
    ]
    ℓweightsₙ = [trace.diagnostics[iter].ℓweightsₙ for iter in eachindex(trace.diagnostics)]
    if isa(smc, _SMC2)
        #Get data and skip latent prediction
        pred = [getfield.(predictionᵛ[iter], 2) for iter in eachindex(predictionᵛ)]
        if length(pred[1][1]) > 1
            return println(
                "Log score only supported for univariate data, use escore instead"
            )
        end
        ℓscoreʷᵉⁱᵍʰᵗᵉᵈ, ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ = compute_ℓscore(
            pred, ℓweightsₙ, data; burnin=burnin
        )
    elseif isa(smc, _IBIS)
        if length(predictionᵛ[1][1]) > 1
            return println(
                "Log score only supported for univariate data, use escore instead"
            )
        end
        ℓscoreʷᵉⁱᵍʰᵗᵉᵈ, ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ = compute_ℓscore(
            predictionᵛ, ℓweightsₙ, data; burnin=burnin
        )
    else
        return println("No escore available for given smc.kernel sampler")
    end

    if (sum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ .== Inf) + sum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ .== Inf)) > 1 && !delete∞
        return println(
            "Non-finite values in log score calculations. If you still want to use log score, use delete∞ = true to replace infinite values with maximal losses",
        )
        #return nothing
    end

    if delete∞
        println(
            "Number of ∞ logpdf evaluations in ℓscoreʷᵉⁱᵍʰᵗᵉᵈ: ",
            sum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ .== Inf),
        )
        println(
            "Number of ∞ logpdf evaluations in ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ: ",
            sum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ .== Inf),
        )
        loss_w = maximum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ[ℓscoreʷᵉⁱᵍʰᵗᵉᵈ .!= Inf])
        loss_uw = maximum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ[ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ .!= Inf])
        println(
            "Replaced non-finite values with max losses, ℓscoreʷ: ",
            loss_w,
            ", ℓscoreʷ",
            loss_uw,
        )
        ℓscoreʷᵉⁱᵍʰᵗᵉᵈ[ℓscoreʷᵉⁱᵍʰᵗᵉᵈ .== Inf] .= loss_w
        ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ[ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ .== Inf] .= loss_uw
    end
    plot_ℓscore = plot(;
        layout=(3, 1),
        size=plot_default_size,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot unnormalized weights
    Plots.scatter!(
        ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ;
        legend=false,
        ylabel="unweighted ℓscore",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=1,
    )
    ## Plot normalized weights
    Plots.scatter!(
        ℓscoreʷᵉⁱᵍʰᵗᵉᵈ;
        legend=false,
        ylabel="weighted ℓscore",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=2,
    )
    ## Plot cumulative versions against unnormalized version
    Plots.plot!(
        cumsum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ);
        label=string("Unweighted ℓscore: ", round(sum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ); digits=2)),
        ylabel="Cumulative ℓscore",
        color="gray",
        legend=:topleft,
        subplot=3,
    )
    Plots.plot!(
        cumsum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ);
        label=string(
            "ℓScores, ", "Weighted ℓscore: ", round(sum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ); digits=2)
        ),
        ylabel="Cumulative ℓscore",
        color="black",
        legend=:topleft,
        subplot=3,
    )
    #xaxis!(string("ℓScores, ", "Weighted ℓscore: ", round(sum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ), digits=2), ", Unweighted ℓscore: ", round(sum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ), digits=2) ), subplot=3)
    xaxis!("Iteration"; subplot=3)
    ## Report ℓscore and return plot
    println("Unweighted ℓscore: ", sum(ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ))
    println("Weighted ℓscore: ", sum(ℓscoreʷᵉⁱᵍʰᵗᵉᵈ))
    display(plot_ℓscore)
    return ℓscoreʷᵉⁱᵍʰᵗᵉᵈ, ℓscoreᵘⁿʷᵉⁱᵍʰᵗᵉᵈ
end

################################################################################
#export
export compute_ℓscore, plot_logscore
