################################################################################
"""
$(SIGNATURES)
Compute escore of trace given data.

# Examples
```julia
```

"""
function compute_escore(prediction::Vector{R}, dataₜ::R) where {R<:Real}
    escore = mean(abs(prediction[iter] - dataₜ) for iter in eachindex(prediction))
    escore -=
        1 / (2 * size(prediction, 1)^2) * sum(
            sum(abs(prediction[main] - prediction[minor])) for main in eachindex(prediction)
            for minor in eachindex(prediction)
        )
    return escore
end
function compute_escore(prediction::Vector{Vector{R}}, dataₜ::T) where {R,T}
    #!NOTE: Euclidean Norm in >1D case, Square root of the sum of the squares
    escore = mean(
        sqrt(sum((prediction[iter] .- dataₜ) .^ 2)) for iter in eachindex(prediction)
    )
    escore -=
        1 / (2 * size(prediction, 1)^2) * sum(
            sum(sqrt(sum((prediction[main] .- prediction[minor]) .^ 2))) for
            main in eachindex(prediction) for minor in eachindex(prediction)
        )
    return escore
end

function compute_escore(predictionᵛ::Vector{Vector{T}}, data::D; burnin=0) where {T,D}
    ## Initiate container
    dataconfig = BaytesCore.ArrayConfig(data)
    #!NOTE: +2 because we start with ( size(latent, 1) - length(predictionᵛ) ) for warmup, then predict for +2 at +1 after warmup
    start_date = size(data, 1) - length(predictionᵛ) + burnin + 2
    #!NOTE: can only go up to T-1, as predictions at endpoint are for T+1, which we do not have.
    escore = zeros(Float64, (length(predictionᵛ) - burnin - 1))
    ## Loop through predictions
    #!NOTE: can only go up to T-1, as predictions at endpoint are for T+1, which we do not have.
    for iter in Base.OneTo((length(predictionᵛ) - burnin - 1))
        ## 1 Forecast to t+1
        predictionₜ₊₁ = predictionᵛ[iter + burnin]
        ## 3 Calculate escores
        #escore[iter]     = compute_escore( predictionₜ₊₁, data[start_date + iter - 1] )
        escore[iter] = compute_escore(
            predictionₜ₊₁, BaytesFilters.grab(data, (start_date + iter - 1), dataconfig)
        )
    end
    ## Return log scores
    return escore
end

################################################################################
"""
$(SIGNATURES)
Plot escore of trace given data.

# Examples
```julia
```

"""
function plot_escore(
    trace::Trace,
    smc,
    data::D;
    burnin=0,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {D}
    ## Get predictions
    predictions = [
        trace.diagnostics[chain].base.prediction for chain in eachindex(trace.diagnostics)
    ]
    if isa(smc, _SMC2)
        #Get data and skip latent prediction
        escore = compute_escore(
            [getfield.(predictions[iter], 2) for iter in eachindex(predictions)],
            data;
            burnin=burnin,
        )
    elseif isa(smc, _IBIS)
        escore = compute_escore(predictions, data; burnin=burnin)
    else
        return println("No escore available for given smc.kernel sampler")
    end
    plot_score = plot(;
        layout=(2, 1),
        size=plot_default_size,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot escore at each iteration
    Plots.scatter!(
        escore;
        legend=false,
        ylabel="Energy-score",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=1,
    )
    ## Plot cumulative e-score
    Plots.plot!(
        cumsum(escore);
        label=string(
            "Cumulative Energy-Score: ",
            round(sum(escore); digits=2),
            " for ",
            size(escore, 1),
            " number of particles. Burnin set to ",
            burnin,
            ".",
        ),
        legend=:topleft,
        ylabel="Cumulative Energy-score",
        xlabel="Time",
        color="gray",
        subplot=2,
    )
    ## Report ℓscore and return escore
    println(
        "Cumulative Energy-Score: ",
        round(sum(escore); digits=2),
        " for ",
        size(escore, 1),
        " number of particles. Burnin set to ",
        burnin,
        ".",
    )
    display(plot_score)
    return escore
end

################################################################################
#export
export compute_escore, plot_escore
