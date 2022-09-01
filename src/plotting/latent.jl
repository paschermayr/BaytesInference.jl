################################################################################
"""
$(SIGNATURES)
Plot latent state trajectory of trace. Only supported for state space models with HMM or HSMM transitions.

# Examples
```julia
```

"""
function plotLatent(trace, algorithm)
    return println("plotLatent not available for current trace and algorithm combination")
end

function plotLatent(
    trace::Trace{C,A,B},
    tagged::Tagged;
    chain=Base.OneTo(length(trace.val)),                              # Choose which chain to plot, by default all chains
    #!FIXME: Come up with better plan to check dimensionality of latent variable, maybe put into pf.tune?
    data=false,                                                       # If "true", plots data below latent states
    latent=false,                                                     # If "true" latent data is assigned, will be plotted alongside estimates
    burnin=0,                                                         # Removes first burnin trajectories from posterior mean calculation
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B} #<:Union{MCMCDiagnostics, PMCMCDiagnostics }}
    val = get_vals(trace, tagged, burnin)
    NTrajectories = length(val[1][1][1][1])
    ## Plot default
    plot_latent = plot(;
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
    ## Plot parameter for each chain
    for Nchain in chain
        ## Calculate Posterior mean in current chain
        if NTrajectories > 1
            θᵗᵉᵐᵖ = [
                vec(
                    mean(
                        getfield.(reduce(hcat, first.(values.(val[Nchain]))), state); dims=2
                    ),
                ) for state in Base.OneTo(NTrajectories)
            ]
        else
            θᵗᵉᵐᵖ = [vec(mean(reduce(hcat, first.(values.(val[Nchain]))); dims=2))]
        end
        for state in Base.OneTo(NTrajectories)
            ## Plot all trajectories ~ too much information
            ## Plot posterior mean
            plot!(
                θᵗᵉᵐᵖ[state];
                lw=1.0,
                palette=Plots.palette(param_color, max(2, NTrajectories)), #color="gold4",
                ylabel=string("latent ", state),
                label=false,
                subplot=state,
            )
        end
    end
    xaxis!("Posterior mean of latent variable for each chain"; subplot=NTrajectories)
    ## Latent - Plot "true" latent data if given
    if latent != false
        if NTrajectories > 1
            latentᵗᵉᵐᵖ = [getfield.(latent, state) for state in Base.OneTo(NTrajectories)]
        else
            latentᵗᵉᵐᵖ = [latent]
        end
        for state in Base.OneTo(NTrajectories)
            plot!(
                latentᵗᵉᵐᵖ[state];
                lw=1.0,
                linestyle=:dot,
                color="black",
                label="Trajectory",
                subplot=state,
            )
        end
    end
    ## Observations
    if data != false
        plot_obs = plot(
            data;
            ylabel="data",
            label=false,
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
        ##Calculate Posterior mean of latent variable for each chain
        if NTrajectories > 1
            θᵐᵉᵃⁿ = [
                [
                    vec(
                        mean(
                            getfield.(reduce(hcat, first.(values.(val[Nchain]))), state);
                            dims=2,
                        ),
                    ) for state in Base.OneTo(NTrajectories)
                ] for Nchain in chain
            ]
            ## Calculate Posterior mean accross chains
            θᵐᵉᵃⁿ⁻ = vec(
                mean(
                    reduce(
                        hcat, [θᵐᵉᵃⁿ[Nchain][1] for Nchain in Base.OneTo(length(chain))]
                    );
                    dims=2,
                ),
            )
        else
            ## Calculate Posterior mean accross chains
            θᵐᵉᵃⁿ = [
                vec(mean(reduce(hcat, first.(values.(val[Nchain]))); dims=2)) for
                Nchain in chain
            ]
            θᵐᵉᵃⁿ⁻ = vec(
                mean(
                    reduce(hcat, [θᵐᵉᵃⁿ[Nchain] for Nchain in Base.OneTo(length(chain))]);
                    dims=2,
                ),
            )
        end
        states_uniform = (θᵐᵉᵃⁿ⁻ .- minimum(θᵐᵉᵃⁿ⁻)) ./ (maximum(θᵐᵉᵃⁿ⁻) - minimum(θᵐᵉᵃⁿ⁻))
        states_obs_scaled = states_uniform .* (maximum(data) - minimum(data)) .+ mean(data)
        states_obs_scaled = states_obs_scaled .- mean(states_obs_scaled) .+ mean(data)
        plot!(
            states_obs_scaled;
            lw=1.0,
            color="gold4",
            label="Rescaled State Posterior Mean of all chains",
        )
        xaxis!("time")#, subplot=2)
    end
    ## Return Plot
    if data == false
        return plot_latent
    else
        return plot(plot_latent, plot_obs; layout=grid(2, 1; heights=[0.7, 0.3]))
    end
end

################################################################################
function plotLatent(
    trace::Trace{C,A,B},
    tagged::Tagged;
    _xaxis = false,
    data=false,                                                           # If "true", plots data below latent states
    latent=false,                                                         # If "true" latent data is assigned, will be plotted alongside estimates
    burnin=0,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B<:SMCDiagnostics}
    ## Unpack Tagged
    val = get_vals(trace, tagged, burnin)
    NTrajectories = length(val[1][1][1][1])
    NIterations = length(val[begin]) #length(trace.val[begin])
    Nparticles = length(trace.val)

    xiter = isa(_xaxis, Bool) ? collect(1:length(val[begin]) ) : _xaxis
    ## Plot default
    plot_latent = plot(;
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
    ## Plot parameter for each chain
    for Nchain in Base.OneTo(NIterations) #(NIterations - burnin))
        ## Calculate Posterior mean in current time step
        if NTrajectories > 1
            θᵗᵉᵐᵖ = [
                vec(
                    mean(
                        getfield.(
                            reduce(
                                hcat,
                                first.(
                                    values.([
                                        (val[particle][Nchain]) for
                                        particle in Base.OneTo(Nparticles)
                                    ]),
                                ),
                            ),
                            state,
                        );
                        dims=2,
                    ),
                ) for state in Base.OneTo(NTrajectories)
            ]
        else
            θᵗᵉᵐᵖ = [
                vec(
                    mean(
                        reduce(
                            hcat,
                            first.(
                                values.([
                                    val[Nparticles][Nchain] for
                                    particle in Base.OneTo(Nparticles)
                                ]),
                            ),
                        );
                        dims=2,
                    ),
                ),
            ]
        end
        for state in Base.OneTo(NTrajectories)
            ## Plot all trajectories ~ too much information
            ## Plot posterior mean
            _xiter = isa(_xaxis, Bool) ? (1:length(θᵗᵉᵐᵖ[state])) : @view(xiter[1:length(θᵗᵉᵐᵖ[state])])
            plot!(
                _xiter,
                θᵗᵉᵐᵖ[state];
                lw=1.0,
#                color=Plots.palette(param_color, (NIterations + 1))[Nchain],
                ylabel=string("latent ", state),
                label=false,
                subplot=state,
            )
        end
    end
    xaxis!(
        "Filtered posterior mean of state trajectory at each time step";
        subplot=NTrajectories,
    )
    ## Latent - Plot "true" latent data if given
    if latent != false
        if NTrajectories > 1
            latentᵗᵉᵐᵖ = [getfield.(latent, state) for state in Base.OneTo(NTrajectories)]
        else
            latentᵗᵉᵐᵖ = [latent]
        end
        for state in Base.OneTo(NTrajectories)
            _xiter = isa(_xaxis, Bool) ? (1:length(latentᵗᵉᵐᵖ[state])) : @view(xiter[1:length(latentᵗᵉᵐᵖ[state])])
            plot!(
                _xiter,
                latentᵗᵉᵐᵖ[state];
                lw=1.0,
                #linestyle=:dot,
                color="black",
                label="Trajectory",
                subplot=state,
            )
        end
    end
    ## Observations
    if data != false
        _xiter = isa(_xaxis, Bool) ? (1:length(data)) : @view(xiter[1:length(data)])
        plot_obs = plot(
            _xiter,
            data;
            ylabel="data",
            label=false,
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
        ## Calculate Posterior mean of latent variable for each chain
        if NTrajectories > 1
            θᵐᵉᵃⁿ = vec(
                mean(
                    getfield.(
                        reduce(
                            hcat,
                            first.(
                                values.([
                                    (val[particle][end]) for
                                    particle in Base.OneTo(Nparticles)
                                ]),
                            ),
                        ),
                        1,
                    );
                    dims=2,
                ),
            )
        else
            θᵐᵉᵃⁿ = vec(
                mean(
                    reduce(
                        hcat,
                        first.(
                            values.([
                                val[particle][end] for particle in Base.OneTo(Nparticles)
                            ]),
                        ),
                    );
                    dims=2,
                ),
            )
        end
        ## Calculate Posterior mean accross chains
        states_uniform = (θᵐᵉᵃⁿ .- minimum(θᵐᵉᵃⁿ)) ./ (maximum(θᵐᵉᵃⁿ) - minimum(θᵐᵉᵃⁿ))
        states_obs_scaled = states_uniform .* (maximum(data) - minimum(data)) #.+ mean( data )
        states_obs_scaled = states_obs_scaled .- mean(states_obs_scaled) .+ mean(data)
        _xiter = isa(_xaxis, Bool) ? (1:length(states_obs_scaled)) : @view(xiter[1:length(states_obs_scaled)])
        plot!(
            _xiter,
            states_obs_scaled;
            lw=1.0,
            color="gold4",
            label="Rescaled State Posterior Mean at final iteration",
        )
        xaxis!("time")#, subplot=2)
    end
    ## Return Plot
    if data == false
        return plot_latent
    else
        return plot(plot_latent, plot_obs; layout=grid(2, 1; heights=[0.7, 0.3]))
    end
end

################################################################################
#export
export plotLatent
