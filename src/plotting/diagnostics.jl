################################################################################
"""
$(SIGNATURES)
Plot samples of trace on individual subplots.

# Examples
```julia
```

"""
################################################################################
# only implemented for univariate data and within SMC so far
function plotDiagnostics(
    trace::Trace{C,A,B},
    transform::TraceTransform;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B}
    return plotDiagnostics(
        trace.diagnostics,
        transform;
        plotsize=plotsize,
        param_color=param_color,
        fontsize=fontsize,
        axissize=axissize,
    )
end
function plotDiagnostics(
    diagnostics::Vector{M},
    transform::TraceTransform;
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {M<:SMCDiagnostics}

    ## Assign Burnin from Transforminfo
    burnin = transform.burnin

    ## Subset diagnostics
    diagnosticsᵛ = diagnostics[1+burnin:end]

    ## Compute rejuvenation time points
    rejuvenated = [diagnosticsᵛ[iter].resampled for iter in eachindex(diagnosticsᵛ)]
    rejuvenated_idx = findall(rejuvenated .== true)
    Niterations = length(diagnosticsᵛ)


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
#        xiter,
        [diagnosticsᵛ[iter].base.temperature for iter in eachindex(diagnosticsᵛ)]; #; subplot=7);
    #    xlims=(1,Niterations),
#        xlims=(xiter[1],xiter[end]),
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

################################################################################
#export
export plotDiagnostics
