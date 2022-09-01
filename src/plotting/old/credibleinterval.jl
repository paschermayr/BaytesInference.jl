################################################################################
"""
$(SIGNATURES)
Plot credible interval for samples of trace.

# Examples
```julia
```

"""
function plotCI(trace, algorithm)
    return println("plotCI not available for current trace and algorithm combination")
end
function plotCI(
    trace::Trace{A,Vector{B}},
    tagged::Tagged;
    chains=Base.OneTo(length(trace.val)),  # Choose which chain to plot, by default all chains
    model=false,                          # If model <: AbstractModel given, plots true parameter
    burnin=0,                             # Removes first burnin parameter from posterior calculations
    CIRegion=[0.025, 0.975],
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {A,B<:Union{MCMCDiagnostics,PMCMCDiagnostics}}
    ## Plot default
    plot_CI = plot(;
        layout=(length(tagged.parameter), 1),
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        size=plot_default_size,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    _lengths = tagged.info.unflatten.unflatten.lengths
    _sz = tagged.info.unflatten.unflatten.sz
    _range = [(_sz[i] - _lengths[i] + 1):_sz[i] for i in eachindex(_sz)]
    ## Plot parameter for each chain
    for Nchain in chains
        ## Calculate Posterior mean and CI
        θ = get_chain(trace.val[Nchain], tagged)
        post_mean = vec(mean(θ[:, (begin + burnin):end]; dims=2))
        CI_temp = collect(
            zip(
                [quantile(θ[i, (begin + burnin):end], CIRegion[1]) for i in 1:size(θ, 1)],
                [quantile(θ[i, (begin + burnin):end], CIRegion[2]) for i in 1:size(θ, 1)],
            ),
        )
        CI = [abs.(CI_temp[iter] .- post_mean[iter]) for iter in eachindex(post_mean)]
        ## Loop through parameter
        for (iter, sym) in enumerate(keys(tagged.parameter))
            ## View to current parameter in total chain via info.positionᵛ, then plot then for each state via info.seq
            post_meanᵗᵉᵐᵖ = @view(post_mean[_range[iter]])
            CIᵗᵉᵐᵖ = @view(CI[_range[iter]])
            Plots.scatter!(
                post_meanᵗᵉᵐᵖ;#[tagged.info[sym].seq[state]]',
                yerror=CIᵗᵉᵐᵖ,#[tagged.info[sym].seq[state]],
                ylabel=sym,
                label="Posterior Mean, 95%CI",
                palette=Plots.palette(param_color, _lengths[iter]),
                subplot=iter,
            )
        end
    end
    xaxis!("Parameter for each chain"; subplot=length(tagged.parameter))
    ## Plot true parameter if model given
    if typeof(model) <: ModelWrapper
        θᵗᵉᵐᵖ = ModelWrappers.flatten(model, tagged)
        for (iter, sym) in enumerate(keys(tagged.parameter))
            θᵛ = θᵗᵉᵐᵖ[_range[iter]]
            Plots.scatter!(
                θᵛ;
                ylabel=sym,
                linestyle=:dot,
                color="gold4",
                shape=:x, #color = Plots.palette(param_color, (Nstates+1) )[state],
                subplot=iter,
            )
        end
    end
    ## Return Plot
    return plot_CI
end

################################################################################
#export
export plotCI
