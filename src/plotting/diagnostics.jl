################################################################################
"""
$(SIGNATURES)
Plot diagnostics of trace.

# Examples
```julia
```

"""
function plotDiagnostics(diagnostics)
    return println("plotDiagnostics not available for diagnostics choice")
end

################################################################################
function plotDiagnostics(
    diagnostics::Vector{Vector{DiagnosticsNUTS}};
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
)
    depth = [
        [diagnostics[chain][iter].depth for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    accept_rate = [
        [
            diagnostics[chain][iter].acceptance_rate for
            iter in eachindex(diagnostics[chain])
        ] for chain in eachindex(diagnostics)
    ]
    ϵ = [
        [diagnostics[chain][iter].ϵ for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    steps = [
        [diagnostics[chain][iter].steps for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    integration_time = [ϵ[chain] .* steps[chain] for chain in eachindex(ϵ)]
    ## Plot default
    plot_diagnostics = plot(;
        layout=(4, 1),
        size=plotsize,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    #plot treedepth
    Plots.plot!(
        reduce(hcat, depth);
        ylabel="Tree depth", #, number of steps = 2^depth", #xlabel = "Iteration for each chain",
        subplot=1,
    )
    #plot average accept rate
    Plots.plot!(
        reduce(hcat, accept_rate);
        ylabel="Avg accept rate",# over whole trajectory", #xlabel = "Iteration for each chain",
        subplot=2,
    )
    #Plot Stepsize
    Plots.plot!(
        log.(reduce(hcat, ϵ));
        ylabel="ℓStepsize", #xlabel = "Iteration for each chain",
        subplot=3,
    )
    #Plot integration time
    Plots.plot!(
        reduce(hcat, integration_time);
        ylabel="Integration time",
        xlabel="Iteration for each chain",
        subplot=4,
    )
    return plot_diagnostics
end
function plotDiagnostics(
    diagnostics::Vector{Vector{DiagnosticsHMC}};
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
)
    ϵ = [
        [diagnostics[chain][iter].ϵ for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    steps = [
        [diagnostics[chain][iter].steps for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    integration_time = [ϵ[chain] .* steps[chain] for chain in eachindex(ϵ)]
    ## Plot default
    plot_diagnostics = plot(;
        layout=(2, 1),
        size=plotsize,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    #Plot Stepsize
    Plots.plot!(
        log.(reduce(hcat, ϵ));
        ylabel="ℓStepsize", #xlabel = "Iteration for each chain",
        subplot=1,
    )
    #Plot integration time
    Plots.plot!(
        reduce(hcat, integration_time);
        ylabel="Integration time",
        xlabel="Iteration for each chain",
        subplot=2,
    )

    return plot_diagnostics
end
function plotDiagnostics(
    diagnostics::Vector{Vector{D}};
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {D<:MCMCKernelDiagnostics}
    ϵ = [
        [diagnostics[chain][iter].ϵ for iter in eachindex(diagnostics[chain])] for
        chain in eachindex(diagnostics)
    ]
    ## Plot default
    plot_diagnostics = plot(;
        layout=(1, 1),
        size=plotsize,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    #Plot Stepsize
    Plots.plot!(
        log.(reduce(hcat, ϵ));
        ylabel="ℓStepsize", #xlabel = "Iteration for each chain",
        subplot=1,
    )
    return plot_diagnostics
end

function plotDiagnostics(
    diagnosticsᵛ::Vector{Vector{M}},
    mcmc::MCMC;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:MCMCDiagnostics}
    ## Plot default
    plot_diagnostics = plot(;
        layout=(2, 1),
        size=plotsize,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot logposterior
    Plots.plot!(
        reduce(
            hcat,
            [
                [diagnosticsᵛ[chain][iter].base.ℓobjective for iter in eachindex(diagnosticsᵛ[chain])]
                for chain in eachindex(diagnosticsᵛ)
            ],
        );
        ylabel="ℓobjective",
        subplot=1,
    )
    ## Plot acceptancerate
    accept = [
        [
            diagnosticsᵛ[chain][iter].accept.accepted for
            iter in eachindex(diagnosticsᵛ[chain])
        ] for chain in eachindex(diagnosticsᵛ)
    ]
    accept_rate = [
        cumsum(accept[iter]) ./ vec(1:length(accept[iter])) for iter in eachindex(accept)
    ]
    Plots.plot!(
        reduce(hcat, accept_rate);
        ylabel="Acceptance rate",
        xlabel="Iteration for each chain",
        subplot=2,
    )
    ## Plot sampler specific diagnostics
    plot_sampler_diagnostics = plotDiagnostics([
        [diagnosticsᵛ[chain][iter].kernel for iter in eachindex(diagnosticsᵛ[chain])] for
        chain in eachindex(diagnosticsᵛ)
    ])
    ## Return main plot
    return plot(
        plot_diagnostics,
        plot_sampler_diagnostics;
        legend=false,
        layout=grid(1, 2; heights=[1.0, 1.0]),
    )
end

################################################################################
function plotDiagnostics(
    diagnosticsᵛ::Vector{Vector{M}},
    pf::ParticleFilter;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:ParticleFilterDiagnostics}
    ## Plot default
    plot_diagnostics = plot(;
        layout=(2, 1),
        size=plotsize,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot likelihood
    return Plots.plot(
        reduce(
            hcat,
            [
                [diagnosticsᵛ[chain][iter].base.ℓobjective for iter in eachindex(diagnosticsᵛ[chain])]
                for chain in eachindex(diagnosticsᵛ)
            ],
        );
        ylabel="ℓlikelihood",
        xlabel="Iteration for each chain",
        size=plot_default_size,
        legend=false,
    )
    #    subplot = 1)
    ## Plot ESS
    #=
    ESS = [ vec( mean(reduce(hcat, diagnostics.ESS[iter]), dims=2) ) for iter in eachindex(diagnostics.ESS) ]
    Plots.plot!(reduce(hcat, ESS),
                ylabel="Average ESS", xlabel = "Iteration for each chain",
                subplot = 2)
    ## Return main plot
    return plot_diagnostics
    =#
end
#=
function plotDiagnostics(
    diagnosticsᵛ::Vector{Vector{M}},
    pmcmc::PMCMC;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:PMCMCDiagnostics}
    mcmc_plot = plotDiagnostics(
        _get_sampler_diagnostics(diagnosticsᵛ, pmcmc),
        pmcmc.kernel.mcmc;
        plotsize=plotsize,
        param_color=param_color,
        fontsize=fontsize,
        axissize=axissize,
    )
    pf_plot = plotDiagnostics(
        _get_sampler_diagnostics(diagnosticsᵛ, pmcmc.kernel.pf),
        pmcmc.kernel.pf;
        plotsize=plotsize,
        param_color=param_color,
        fontsize=fontsize,
        axissize=axissize,
    )
    return plot(mcmc_plot, pf_plot; legend=false, layout=grid(2, 1; heights=[0.70, 0.30]))
end
=#
################################################################################

function plotDiagnostics(
    diagnosticsᵛ::Vector{M},
    smc;
    _xaxis = false,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:SMCDiagnostics}
    ## Compute rejuvenation time points
    rejuvenated = [diagnosticsᵛ[iter].resampled for iter in eachindex(diagnosticsᵛ)]
    rejuvenated_idx = findall(rejuvenated .== true)
    Niterations = length(diagnosticsᵛ)

    xiter = isa(_xaxis, Bool) ? collect(1:Niterations) : _xaxis[end-Niterations+1:end]


    ## Plot default
    plot_diagnostics = plot(;
        layout=(7, 1),
        size=plotsize,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot ℓlikelihood
    Plots.scatter!(
        reduce(hcat, [diagnosticsᵛ[iter].ℓweights for iter in eachindex(diagnosticsᵛ)])';
    #    xaxis=false,
        xticks=false,
        xlims=(1,Niterations),
        ylabel="ℓℒ estimates",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=1,
    )

    ## Plot ℓlikelihood variance
    Plots.plot!(
        std.([diagnosticsᵛ[iter].ℓweights for iter in eachindex(diagnosticsᵛ)]) .^ 2;
    #    xaxis=false,
        xticks=false,
        xlims=(1,Niterations),
        ylabel="ℓℒ estimate variance",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=2,
    )
    ## Plot theta correlation
    Plots.scatter!(
        rejuvenated_idx,
        reduce(
            hcat,
            [diagnosticsᵛ[rejuvenated][iter].ρ for iter in Base.OneTo(sum(rejuvenated))],
        )';
    #    xaxis=false,
        xticks=false,
        xlims=(1,Niterations),
        ylabel="Rej. ρ for θ",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=3,
    )
    #Plot average of theta correlation at each step
    Plots.scatter!(
        rejuvenated_idx,
        mean.([diagnosticsᵛ[rejuvenated][iter].ρ for iter in Base.OneTo(sum(rejuvenated))]);
    #    xaxis=false,
        xticks=false,
        #ylabel="Rej. ρ for θ",
        xlims=(1,Niterations),
        color="red",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:x,
        markerstrokecolor="grey",
        markersize=2,
        subplot=3,
    )
    ## Plot number of jittersteps
    if length(rejuvenated_idx) > 1
        Plots.scatter!(
            rejuvenated_idx,
            reduce(
                hcat,
                [
                    diagnosticsᵛ[rejuvenated][iter].jittersteps for
                    iter in Base.OneTo(sum(rejuvenated))
                ],
            )';
        #    xaxis=false,
            xticks=false,
            xlims=(1,Niterations),
            ylabel="Rej. steps",
            color="black",
            markerstrokewidth=0.0,
            alpha=1.0,
            shape=:o,
            markerstrokecolor="grey",
            markersize=2,
            subplot=4,
        )

    #Plot average of theta correlation at each step
    Plots.scatter!(
        rejuvenated_idx,
        mean.([
            diagnosticsᵛ[rejuvenated][iter].jittersteps for
            iter in Base.OneTo(sum(rejuvenated))
        ]);
    #    xaxis=false,
        xticks=false,
        #ylabel="Rej. ρ for θ",
        xlims=(1,Niterations),
        color="red",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:x,
        markerstrokecolor="grey",
        markersize=2,
        subplot=4,
    )
end
    ## Plot normalized ℓweights
    Plots.scatter!(
        exp.(
            reduce(hcat, [diagnosticsᵛ[iter].ℓweightsₙ for iter in eachindex(diagnosticsᵛ)])
        )';
    #    xaxis=false,
        xticks=false,
        xlims=(1,Niterations),
        ylabel="Norm. weights after Rej.",
        color="black",
        markerstrokewidth=0.0,
        alpha=1.0,
        shape=:o,
        markerstrokecolor="grey",
        markersize=2,
        subplot=5,
    )
    ## Plot ESS
    Plots.plot!(
        [diagnosticsᵛ[iter].ESS for iter in eachindex(diagnosticsᵛ)];
    #    xaxis=false,
    #    xticks=false,
        xlims=(1,Niterations),
        ylabel="ESS before Rej.",
        color="black",
        subplot=6,
    )
    ## Plot Temperature
    Plots.plot!(
        xiter,
        [diagnosticsᵛ[iter].base.temperature for iter in eachindex(diagnosticsᵛ)]; #; subplot=7);
    #    xlims=(1,Niterations),
        xlims=(xiter[1],xiter[end]),
        ylabel="Temperature",
        color="black",
        subplot=7,
    )
    ## Plot accepted - do not need anymore is easily visible by other graphs
    #=
        Plots.scatter!(rejuvenated_idx, repeat([1], length(rejuvenated_idx) ),
                    markersize = 2,
                    ylabel="θ Rej. times", xlabel = "Iteration for each particle, Rej. = Particle rejuvenation",
                    subplot = 6)
    =#
    xaxis!(
        "Iteration for each particle, Rej. = Particle rejuvenation";
        subplot=size(plot_diagnostics, 1),
    )
    ## Return main plot
    return plot_diagnostics
end

function plotDiagnostics(trace::Trace, algorithm; kwargs...)
    return plotDiagnostics(trace.diagnostics, algorithm; kwargs...)
end

################################################################################
#export
export plotDiagnostics
