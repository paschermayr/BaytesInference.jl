################################################################################
"""
$(SIGNATURES)
Plot samples of trace on individual subplots.

# Examples
```julia
```

"""
function plotCredibleInterval(
    trace::Trace{C,A,B},
    transform::TraceTransform;
    model=false,                            # If model <: AbstractModel given, plots true parameter
    dates = false,
    paramnames = transform.paramnames,          # Get Model names
    CIRegion = [0.025, 0.975], 
    layout = length_constrained(transform.tagged),
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B}
    ## Plot default
    plot_chains = plot(;
        layout= layout, #(length_constrained(transform.tagged), 1),
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        size=plotsize,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Get Trace Values
    _vals = trace_to_3DArray(trace, transform)
    _Niter = size(_vals, 1)
    _Nchains = size(_vals, 2)
    _Nparam = size(_vals, 3)
    burnin = transform.burnin
    maxiterations = transform.maxiterations

    # Obtain X scale and check for matching indices
    xiter = dates == false ? collect(1:maxiterations)[(burnin+1):end] : dates[end-_Niter+1:end]
    @argcheck length(paramnames) == _Nparam "Number of Parameter subset and Parameter Names do not match"
    @argcheck length(xiter) == _Niter "Dates has different index than chain"

    ## if model there, obtain true parameter
    if model != false
        θ_true = flatten(model, transform.tagged)
        @argcheck length(θ_true) == size(_vals, 3)
    end
    ## Plot
    for iter in eachindex(paramnames)
        # Obtain All samples for current parameter
        θₜₑₘₚ = view(_vals, :, :, iter)
        # Compute Posterior Mean and Credible Interval
        post_mean = mean(θₜₑₘₚ, dims=2)
        CI = Vector{Tuple{Float64,Float64}}(undef, _Niter)
        for idx in eachindex(CI)
            θᵢ = θₜₑₘₚ[idx, :] 
            CI[idx] = (
                abs.(quantile(θᵢ, CIRegion[1]) - post_mean[idx]),
                abs.(quantile(θᵢ, CIRegion[2]) - post_mean[idx]),
            )
        end
        ## Plot Credible Interval
        Plots.plot!(xiter, post_mean;## Posterior Mean
            shape=:o,
            markerstrokewidth=0.1,
            markersize=2,
            ylabel= paramnames[iter],
            label=false, #string("Posterior Mean, CI: (", CIRegion[1],", ", CIRegion[2], ")"),
            color="black", #"gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
            ## CI
            #!NOTE: Internally Matrices just viewed as vector of vectors in row major order, so can just unroll Matrix CI Tuples as Tuple of individual vectors
            ribbon=(
                ( getfield.(CI, 1), getfield.(CI, 2) )
            ),
            fillcolor="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
            fillalpha=0.25,
            subplot = iter,
        )
        ## If defined, plot true parameter as vline
        if model != false
            Plots.hline!(
                [θ_true[iter]];
                label = false,
                ylabel = paramnames[iter],
                linestyle=:dot,
                color = "black", #palette = Plots.palette(param_color, _Nchains),
                subplot=iter,
            )
        end
    end
    return plot_chains
end

################################################################################
#export
export plotCredibleInterval
