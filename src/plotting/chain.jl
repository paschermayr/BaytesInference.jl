################################################################################
"""
$(SIGNATURES)
Plot samples of trace.

# Examples
```julia
```

"""
function plotChain(trace, algorithm)
    return println("plotChain not available for current trace and algorithm combination")
end

function plotChain(
    trace::Trace{C,A,B},
    tagged::Tagged; 
    chains=Base.OneTo(length(trace.val)),  # Choose which chain to plot, by default all chains
     model=false,                          # If model <: AbstractModel given, plots true parameter
    burnin=0,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B}
    ## Plot default
    plot_chain = plot(;
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
    ## Get indices for individual parameter
    _lengths = tagged.info.reconstruct.unflatten.strict._unflatten.lngth
    _sz = tagged.info.reconstruct.unflatten.strict._unflatten.sz
    _range = [(_sz[i] - _lengths[i] + 1):_sz[i] for i in eachindex(_sz)]
    ## Plot parameter for each chain - I think I should just flatten everything here?
    θ = get_vals(trace, tagged, burnin) #Need it that way because from flatten_chains I dont get dimensionality
    for Nchain in chains
        for (iter, sym) in enumerate(keys(tagged.parameter))
            reconstruct = ModelWrappers.ReConstructor(subset(θ[Nchain][begin], sym))
            θ_temp = reduce(
                hcat,
                flatten(reconstruct, subset(θ[Nchain][idx], sym)) for idx in eachindex(θ[Nchain])
            )
            ## view to current parameter in total chain via info.positionᵛ, then plot then for each state via info.seq
            Plots.plot!(
                θ_temp';
                ylabel=sym,
                palette=if _lengths[iter] > 1
                    Plots.palette(param_color, _lengths[iter])
                else
                    Plots.palette(param_color, 2)
                end,
                subplot=iter,
            )
        end
    end
    xaxis!("MCMC iteration for each chain"; subplot=length(tagged.parameter))
    ## Plot true parameter if model given
    if typeof(model) <: ModelWrapper
        θᵗᵉᵐᵖ = ModelWrappers.flatten(model, tagged)
        for (iter, sym) in enumerate(keys(tagged.parameter))
            θᵛ = θᵗᵉᵐᵖ[_range[iter]]
            Plots.hline!(
                θᵛ';
                ylabel=sym,
                linestyle=:dot,
                palette=Plots.palette(param_color, _lengths[iter]),
                subplot=iter,
            )
        end
    end
    ## Return Plot
    return plot_chain
end

################################################################################
function plotChain(
    trace::Trace{C,A,B},
    tagged::Tagged;
    _xaxis = false,
    ShowCI=true,
    CIRegion=[0.025, 0.975],
    model=false,                          # If model <: AbstractModel given, plots true parameter
    burnin=0,
    flattenparam=FlattenDefault(),
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B<:SMCDiagnostics}
    ## Get indices for individual parameter
    _lengths = tagged.info.reconstruct.unflatten.strict._unflatten.lngth #tagged.info.unflatten.unflatten.lengths
    _sz = tagged.info.reconstruct.unflatten.strict._unflatten.sz
    _range = [(_sz[i] - _lengths[i] + 1):_sz[i] for i in eachindex(_sz)]

    N = length(trace.val[begin]) - burnin
    post_mean = zeros(Float64, length(tagged), N)
    CI = Matrix{Tuple{Float64,Float64}}(undef, length(tagged), N)

    for iter in Base.OneTo(N)
        θᵗᵉᵐᵖ = flatten_index(trace, tagged, iter + burnin, flattenparam) #get_chain(trace.val[iter+burnin], tagged)
        post_mean[:, iter] = vec(mean(θᵗᵉᵐᵖ; dims=2))
        CI_temp = collect(
            zip(
                [quantile(θᵗᵉᵐᵖ[i, :], CIRegion[1]) for i in 1:size(θᵗᵉᵐᵖ, 1)],
                [quantile(θᵗᵉᵐᵖ[i, :], CIRegion[2]) for i in 1:size(θᵗᵉᵐᵖ, 1)],
            ),
        )
        CI[:, iter] = [abs.(CI_temp[i] .- post_mean[i, iter]) for i in eachindex(CI_temp)]
    end
    xiter = isa(_xaxis, Bool) ? collect(1:size(post_mean, 2) ) : _xaxis[end-size(post_mean, 2)+1:end]

    ## Plot default
    plot_chain = plot(;
        layout=(length(tagged.parameter), 1),
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
    ## Loop through parameter
    for (iter, sym) in enumerate(keys(tagged.parameter))
        ## View to current parameter in total chain via info.positionᵛ, then plot then for each state via info.seq
        post_meanᵗᵉᵐᵖ = @view(post_mean[_range[iter], :])
        CIᵗᵉᵐᵖ = @view(CI[_range[iter], :])
        if ShowCI
            for θdim in Base.OneTo(size(post_meanᵗᵉᵐᵖ, 1)) #Necessary for correct Credible Intervals in Multivariate Case
                Plots.plot!(
                    xiter,
                    post_meanᵗᵉᵐᵖ[θdim, :];
                    ## Posterior Mean
                    shape=:o,
                    markerstrokewidth=0.1,
                    markersize=2,
                    ylabel=sym,
                    label="Posterior Mean, 95%CI",
                    #                palette = Plots.palette(param_color, _lengths[iter]),
                    palette = if _lengths[iter] > 1
                        Plots.palette(param_color, _lengths[iter])
                    else
                        Plots.palette(param_color, 2)
                    end,
                    ## CI
                    #!NOTE: Internally Matrices just viewed as vector of vectors in row major order, so can just unroll Matrix CI Tuples as Tuple of individual vectors
                    ribbon=(
                        getfield.(vec(CIᵗᵉᵐᵖ[θdim, :]), 1),
                        getfield.(vec(CIᵗᵉᵐᵖ[θdim, :]), 2),
                    ),
                    #    fillpalette = Plots.palette(param_color, _lengths[iter]),
                    fillpalette = if _lengths[iter] > 1
                        Plots.palette(param_color, _lengths[iter])
                    else
                        Plots.palette(param_color, 2)
                    end,
                    fillalpha=0.25,
                    #markerstrokewidth=0.1, alpha = 1.0,
                    subplot=iter,
                )
            end
        else
            Plots.plot!(
                xiter,
                post_meanᵗᵉᵐᵖ';
                ## Posterior Mean
                shape=:o,
                markerstrokewidth=0.1,
                markersize=2,
                ylabel=sym,
                label="Posterior Mean, 95%CI",
                palette=Plots.palette(param_color, _lengths[iter]),
                subplot=iter,
            )
        end
        #end
    end
    xaxis!(
        string(
            "MCMC iteration for each chain",
            " - Posterior Mean, CI: (",
            CIRegion[1],
            ", ",
            CIRegion[2],
            ")",
        );
        subplot=length(tagged.parameter),
    )
    ## Plot true parameter if model given
    if typeof(model) <: ModelWrapper
        θᵗᵉᵐᵖ = ModelWrappers.flatten(model, tagged)
        for (iter, sym) in enumerate(keys(tagged.parameter))
            θᵛ = θᵗᵉᵐᵖ[_range[iter]]
            Plots.hline!(
                θᵛ';
                ylabel=sym,
                linestyle=:dot,
                palette=Plots.palette(param_color, _lengths[iter]),
                subplot=iter,
            )
        end
    end
    ## Return Plot
    return plot_chain
end

################################################################################
#export
export plotChain
