################################################################################
"""
$(SIGNATURES)
Plot samples of trace on individual subplots.

# Examples
```julia
```

"""
function plotChains(
    trace::Trace{C,A,B},
    transform::TraceTransform;
    model=false,                            # If model <: AbstractModel given, plots true parameter
    dates = false,
    paramnames = transform.paramnames,          # Get Model names
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
    xiter = dates == false ? collect(1:maxiterations)[(burnin+1):end] : dates[(burnin+1):end]
    @argcheck length(paramnames) == _Nparam "Number of Parameter subset and Parameter Names do not match"
    @argcheck length(xiter) == _Niter "Dates has different index than chain"


    ## if model there, obtain true parameter
    if model != false
        θ_true = flatten(model, transform.tagged)
        @argcheck length(θ_true) == size(_vals, 3)
    end
    ## Plot
    for iter in eachindex(paramnames)
        plot!(xiter, view(_vals, :, :, iter),
            label = false,
            ylabel = paramnames[iter],
            palette = Plots.palette(param_color, _Nchains),
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
export plotChains
