################################################################################
"""
$(SIGNATURES)
Plot posterior predictive samples of trace. Only supported for some algorithms.

# Examples
```julia
```

"""
function plotPosteriorPrediction(trace, algorithm)
    return println(
        "plotPosteriorPrediction not available for current trace and algorithm combination"
    )
end

function plotPosteriorPrediction(
    diagnosticsᵛ::Vector{Vector{M}},
    pf::ParticleFilter;
    chain=Base.OneTo(length(diagnosticsᵛ)),                                  # Choose which chain to plot, by default all chains
    burnin=0,                                                             # Removes first burnin trajectories from posterior mean calculation
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:ParticleFilterDiagnostics}
    ## Get predictions
    predictions = [
        [diagnosticsᵛ[chain][iter].base.prediction for iter in eachindex(diagnosticsᵛ[chain])]
        for chain in eachindex(diagnosticsᵛ)
    ]
    NTrajectories = length(predictions[1][1][1])   # Guess the dimensionality of trajectory
    Ndata = size(predictions[1][1][2], 1)
    ## Plot default
    plot_prediction = plot(;
        layout=(NTrajectories, 1),
        size=plot_default_size,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=:topright,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot state for each chain
    for Nchain in chain
        ## Obtain all states in current chain
        if NTrajectories > 1
            #Obtain states first, then split states
            stateᵗ⁺¹ = [
                getfield.(getfield.(predictions[Nchain][(begin + burnin):end], 1), state)
                for state in Base.OneTo(NTrajectories)
            ]
        else
            stateᵗ⁺¹ = [getfield.(predictions[Nchain][(begin + burnin):end], 1)]
        end
        for state in Base.OneTo(NTrajectories)
            ## Plot all trajectories ~ too much information
            ## Plot posterior mean
            histogram!(
                stateᵗ⁺¹[state];
                alpha=0.5,
                color=Plots.palette(param_color, (length(chain) + 1))[Nchain], #color="gold4",
                ylabel=string("Latent state ", state),
                label=string("Chain ", Nchain),
                subplot=state,
            )
        end
    end
    xaxis!("Posterior predictive distribution for State"; subplot=1)
    if NTrajectories > 1
        xaxis!("Time remaining in predicted state"; subplot=NTrajectories)
    end
    ## Plot data for each chain
    plot_obs = plot(;
        ylabel="data",
        xlabel="Posterior predictive distribution for data given state",
        label=false,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        size=plot_default_size,
        legend=:topright,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    for Nchain in chain
        ## Obtain all data in current chain
        dataᵗ⁺¹ = getfield.(predictions[Nchain][(begin + burnin):end], 2)
        if !isa(dataᵗ⁺¹, Vector{R where R<:Real})
            dataᵗ⁺¹ = reduce(hcat, dataᵗ⁺¹)'
        end
        for datadim in Base.OneTo(Ndata)
            histogram!(
                dataᵗ⁺¹[:, datadim];
                alpha=0.25,
                color=Plots.palette(param_color, (length(Ndata) + 1))[datadim], #color="gold4",
                label=string("Data dim ", datadim, ", Chain ", Nchain),
            )
        end
    end
    ## Return Plot
    return plot(plot_prediction, plot_obs; layout=grid(2, 1; heights=[0.7, 0.3]))
end

################################################################################
function plotPosteriorPrediction(
    diagnosticsᵛ::Vector{Vector{M}},
    mcmc::MCMC;
    chain=Base.OneTo(length(diagnosticsᵛ)),                                  # Choose which chain to plot, by default all chains
    burnin=0,                                                             # Removes first burnin trajectories from posterior mean calculation
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:MCMCDiagnostics}
    ## Get predictions
    predictions = [
        [diagnosticsᵛ[chain][iter].base.prediction for iter in eachindex(diagnosticsᵛ[chain])]
        for chain in eachindex(diagnosticsᵛ)
    ]
    ## Plot data for each chain
    Ndata = size(predictions[1][1], 1)
    plot_obs = plot(;
        ylabel="data",
        xlabel="Posterior predictive distribution",
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        label=false,
        size=plot_default_size,
        legend=:topright,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    for Nchain in chain
        ## Obtain all data in current chain
        dataᵗ⁺¹ = predictions[Nchain][(begin + burnin):end]
        if !isa(dataᵗ⁺¹, Vector{R where R<:Real})
            dataᵗ⁺¹ = reduce(hcat, dataᵗ⁺¹)'
        end
        for datadim in Base.OneTo(Ndata)
            histogram!(
                dataᵗ⁺¹[:, datadim];
                alpha=0.25,
                color=Plots.palette(param_color, (length(Ndata) + 1))[datadim], #color="gold4",
                label=string("Data dim ", datadim, ", Chain ", Nchain),
            )
        end
    end
    ## Return Plot
    return plot_obs
end
################################################################################
#=
function plotPosteriorPrediction(
    diagnosticsᵛ::Vector{Vector{M}},
    pmcmc::PMCMC;
    chain=Base.OneTo(length(diagnosticsᵛ)),                                  # Choose which chain to plot, by default all chains
    burnin=0,                                                             # Removes first burnin trajectories from posterior mean calculation
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:PMCMCDiagnostics}
    pf_plot = plotPosteriorPrediction(
        _get_sampler_diagnostics(diagnosticsᵛ, pmcmc.kernel.pf),
        pmcmc.kernel.pf;
        chain=chain,
        burnin=burnin,
        plotsize=plotsize,
        param_color=param_color,
        fontsize=fontsize,
        axissize=axissize,
    )
    mcmc_plot = plotPosteriorPrediction(
        _get_sampler_diagnostics(diagnosticsᵛ, pmcmc.kernel.mcmc),
        pmcmc.kernel.mcmc;
        chain=chain,
        burnin=burnin,
        plotsize=plotsize,
        param_color=param_color,
        fontsize=fontsize,
        axissize=axissize,
    )
    return plot(
        pf_plot,
        mcmc_plot; #legend=false,
        layout=grid(2, 1; heights=[0.80, 0.20]),
    )
end
=#
################################################################################
function plotPosteriorPrediction(
    diagnosticsᵛ::Vector{M},
    smc;
    _xaxis = false,
    data=false,                                                           # If "true", plots data below latent states
    latent=false,                                                         # If "true" latent data is assigned, will be plotted alongside estimates
    burnin=0,
    ShowCI=true,
    CIRegion=[0.025, 0.975],
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:SMCDiagnostics}

    ################################################################################
    # SMC2
    ################################################################################
    ################################################# LATENT DATA
    if isa(smc, _SMC2)
        ## Grab predictions
        predictions = [diagnosticsᵛ[chain].base.prediction for chain in eachindex(diagnosticsᵛ)]
        NTrajectories = length(predictions[1][1][1])   # Guess the dimensionality of trajectory
        Ndata = size(predictions[1][1][2], 1)
        N = length(predictions) - burnin
        xiter = isa(_xaxis, Bool) ? collect(1:length(predictions) ) : _xaxis[end-length(predictions)+1:end]
        ## Plot default
        plot_latent = plot(;
            #xiter,
            layout=(NTrajectories, 1),
            size=plot_default_size,
            foreground_color_legend = :transparent,
            background_color_legend = :transparent,
            legend=:topleft,
            xguidefontsize=fontsize,
            yguidefontsize=fontsize,
            legendfontsize=fontsize,
            xtickfontsize=axissize,
            ytickfontsize=axissize,
        )
        ## Create Posterior mean and CI for each iteration
        post_mean = zeros(Float64, NTrajectories, N)
        CI = Matrix{Tuple{Float64,Float64}}(undef, NTrajectories, N)
        for Nchain in Base.OneTo(N)
            if NTrajectories > 1
                stateᵗ⁺¹ = [
                    getfield.(getfield.(predictions[Nchain + burnin], 1), state) for
                    state in Base.OneTo(NTrajectories)
                ]
            else
                stateᵗ⁺¹ = [getfield.(predictions[Nchain + burnin], 1)]
            end
            for state in Base.OneTo(NTrajectories)
                stateᵗᵉᵐᵖ = stateᵗ⁺¹[state]
                post_mean[state, Nchain] = mean(stateᵗᵉᵐᵖ)
                CI[state, Nchain] = (
                    post_mean[state, Nchain] - quantile(stateᵗᵉᵐᵖ, CIRegion[1]),
                    quantile(stateᵗᵉᵐᵖ, CIRegion[2]) - post_mean[state, Nchain],
                )
            end
        end
        ## Plot CI
        for state in Base.OneTo(NTrajectories)
            if ShowCI
                Plots.plot!(
                    xiter[end-length(post_mean[state, :])+1:end],
                    post_mean[state, :];
                    ## Posterior Mean
                    shape=:o,
                    markerstrokewidth=0.1,
                    markersize=2,
                    ylabel=string("Latent state ", state),
                    label=false, #string("Posterior Mean, CI: (", CIRegion[1],", ", CIRegion[2], ")"),
                    color="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    ## CI
                    #!NOTE: Internally Matrices just viewed as vector of vectors in row major order, so can just unroll Matrix CI Tuples as Tuple of individual vectors
                    ribbon=(
                        getfield.(vec(CI[state, :]), 1), getfield.(vec(CI[state, :]), 2)
                    ),
                    fillcolor="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    fillalpha=0.25,
                    #markerstrokewidth=0.1, alpha = 1.0,
                    subplot=state,
                )
            else
                Plots.plot!(
                    xiter[end-length(post_mean[state, :])+1:end],
                    post_mean[state, :];
                    ## Posterior Mean
                    shape=:o,
                    markerstrokewidth=0.1,
                    markersize=2,
                    ylabel=string("Latent state ", state),
                    label=false, #string("Posterior Mean, CI: (", CIRegion[1],", ", CIRegion[2], ")"),
                    color="gold4",
                    subplot=state,
                )
            end
        end
        xaxis!("Posterior predictive mean for t+1"; subplot=NTrajectories)
        ## Plot against sample latent data if available
        if latent != false
            #!NOTE: +2 because we start with ( size(latent, 1) - length(predictionᵛ) ) for warmup, then predict for +2 at +1 after warmup
            start_date = size(latent, 1) - length(predictions) + burnin + 2 #Should have 1 data point less than chains above

            if NTrajectories > 1
                latentᵗᵉᵐᵖ = [
                    getfield.(latent, state) for state in Base.OneTo(NTrajectories)
                ]
            else
                latentᵗᵉᵐᵖ = [latent]
            end
            for state in Base.OneTo(NTrajectories)
                Niter = length(latentᵗᵉᵐᵖ[state][start_date:end])
                Plots.plot!(
                    xiter[end-Niter+1:end],
                    latentᵗᵉᵐᵖ[state][start_date:end];
                    linewidth=0.5,
                    label="Sample data",
                    color="black", #linestyle=:dot,
                    subplot=state,
                )
            end
        end
        ################################################# OBSERVED DATA
        ## Create Posterior mean and CI for each iteration
        post_mean = zeros(Float64, Ndata, N)
        CI = Matrix{Tuple{Float64,Float64}}(undef, Ndata, N)
        for Nchain in Base.OneTo(N)
            dataᵗ⁺¹ = getfield.(predictions[Nchain + burnin], 2)
            if !isa(dataᵗ⁺¹, Vector{R} where {R<:Real})
                dataᵗ⁺¹ = reduce(hcat, dataᵗ⁺¹)'
            end
            for datadim in Base.OneTo(Ndata)
                dataᵗᵉᵐᵖ = dataᵗ⁺¹[:, datadim]
                post_mean[datadim, Nchain] = mean(dataᵗᵉᵐᵖ)
                CI[datadim, Nchain] = (
                    abs.(quantile(dataᵗᵉᵐᵖ, CIRegion[1]) - post_mean[datadim, Nchain]),
                    abs.(quantile(dataᵗᵉᵐᵖ, CIRegion[2]) - post_mean[datadim, Nchain]),
                )
            end
        end
        ## Plot data
        plot_data = plot(;
            layout=(Ndata, 1),
            size=plot_default_size,
            foreground_color_legend = :transparent,
            background_color_legend = :transparent,
            xguidefontsize=fontsize,
            yguidefontsize=fontsize,
            legendfontsize=fontsize,
            xtickfontsize=axissize,
            ytickfontsize=axissize,
        )
        ## Plot CI
        for datadim in Base.OneTo(Ndata)
            if ShowCI
                Plots.plot!(
                    xiter[end-length(post_mean[datadim, :])+1:end],
                    post_mean[datadim, :];
                    ## Posterior Mean
                    shape=:o,
                    markerstrokewidth=0.1,
                    markersize=2,
                    ylabel=string("Data dimension ", datadim),
                    label=false, #string("Posterior Mean, CI: (", CIRegion[1],", ", CIRegion[2], ")"),
                    color="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    ## CI
                    #!NOTE: Internally Matrices just viewed as vector of vectors in row major order, so can just unroll Matrix CI Tuples as Tuple of individual vectors
                    ribbon=(
                        getfield.(vec(CI[datadim, :]), 1), getfield.(vec(CI[datadim, :]), 2)
                    ),
                    fillcolor="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    fillalpha=0.25,
                    #markerstrokewidth=0.1, alpha = 1.0,
                    subplot=datadim,
                )
            end
        end
        xaxis!(
            string(
                "Posterior predictive mean for t+1",
                " - Posterior Mean, CI: (",
                CIRegion[1],
                ", ",
                CIRegion[2],
                ")",
            );
            subplot=Ndata,
        )
        ## Plot against sample latent data if available
        if data != false
            #!NOTE: +2 because we start with ( size(latent, 1) - length(predictionᵛ) ) for warmup, then predict for +2 at +1 after warmup
            start_date = size(data, 1) - length(predictions) + burnin + 2 #Should have 1 data point less than chains above
            for datadim in Base.OneTo(Ndata)
                Niter = length(data[start_date:end, datadim])
                Plots.plot!(
                    xiter[end-Niter+1:end],
                    data[start_date:end, datadim];
                    linewidth=0.5,
                    label="Sample data",
                    color="black", #linestyle=:dot,
                    subplot=datadim,
                )
            end
        end
        ## Return Plot
        return plot(plot_latent, plot_data; layout=grid(2, 1; heights=[0.5, 0.5]))
        ################################################################################
        # IBIS
        ################################################################################
    elseif isa(smc, _IBIS)
        ## Grab predictions
        predictions = [diagnosticsᵛ[chain].prediction for chain in eachindex(diagnosticsᵛ)]
        Ndata = size(predictions[1][1], 1)
        N = length(predictions) - burnin
        xiter = isa(_xaxis, Bool) ? collect(1:length(predictions) ) : _xaxis[end-length(predictions)+1:end]

        ################################################# OBSERVED DATA
        ## Create Posterior mean and CI for each iteration
        post_mean = zeros(Float64, Ndata, N)
        CI = Matrix{Tuple{Float64,Float64}}(undef, Ndata, N)
        for Nchain in Base.OneTo(N)
            dataᵗ⁺¹ = predictions[Nchain + burnin]
            if !isa(dataᵗ⁺¹, Vector{R} where {R<:Real})
                dataᵗ⁺¹ = reduce(hcat, dataᵗ⁺¹)'
            end
            for datadim in Base.OneTo(Ndata)
                dataᵗᵉᵐᵖ = dataᵗ⁺¹[:, datadim]
                post_mean[datadim, Nchain] = mean(dataᵗᵉᵐᵖ)
                CI[datadim, Nchain] = (
                    abs.(quantile(dataᵗᵉᵐᵖ, CIRegion[1]) - post_mean[datadim, Nchain]),
                    abs.(quantile(dataᵗᵉᵐᵖ, CIRegion[2]) - post_mean[datadim, Nchain]),
                )
            end
        end
        ## Plot data
        plot_data = plot(;
            layout=(Ndata, 1),
            size=plot_default_size,
            foreground_color_legend = :transparent,
            background_color_legend = :transparent,
            legend=:topleft,
            xguidefontsize=fontsize,
            yguidefontsize=fontsize,
            legendfontsize=fontsize,
            xtickfontsize=axissize,
            ytickfontsize=axissize,
        )
        ## Plot CI
        for datadim in Base.OneTo(Ndata)
            if ShowCI
                Plots.plot!(
                    xiter[end-length(post_mean[datadim, :])+1:end],
                    post_mean[datadim, :];
                    ## Posterior Mean
                    shape=:o,
                    markerstrokewidth=0.1,
                    markersize=2,
                    ylabel=string("Data dimension ", datadim),
                    label=false, #string("Posterior Mean, CI: (", CIRegion[1],", ", CIRegion[2], ")"),
                    color="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    ## CI
                    #!NOTE: Internally Matrices just viewed as vector of vectors in row major order, so can just unroll Matrix CI Tuples as Tuple of individual vectors
                    ribbon=(
                        getfield.(vec(CI[datadim, :]), 1), getfield.(vec(CI[datadim, :]), 2)
                    ),
                    fillcolor="gold4", #Plots.palette(param_color, (NTrajectories+1) )[state],
                    fillalpha=0.25,
                    #markerstrokewidth=0.1, alpha = 1.0,
                    subplot=datadim,
                )
            end
        end
        xaxis!(
            string(
                "Posterior predictive mean for t+1",
                " - Posterior Mean, CI: (",
                CIRegion[1],
                ", ",
                CIRegion[2],
                ")",
            );
            subplot=Ndata,
        )
        ## Plot against sample latent data if available
        if data != false
            #!NOTE: +2 because we start with ( size(latent, 1) - length(predictionᵛ) ) for warmup, then predict for +2 at +1 after warmup
            start_date = size(data, 1) - length(predictions) + burnin + 2 #Should have 1 data point less than chains above
            for datadim in Base.OneTo(Ndata)
                Niter = length(data[start_date:end, datadim])
                Plots.plot!(
                    xiter[end-Niter+1:end],
                    data[start_date:end, datadim];
                    linewidth=0.5,
                    label="Sample data",
                    color="black", #linestyle=:dot,
                    subplot=datadim,
                )
            end
        end
        ## Return Plot
        return plot_data
    else
        println("No plotting function for smc.kernel available")
    end
end

function plotPosteriorPrediction(trace::Trace, algorithm)
    return plotPosteriorPrediction(trace.diagnostics, algorithm)
end

################################################################################
#export
export plotPosteriorPrediction
