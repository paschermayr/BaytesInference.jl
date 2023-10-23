################################################################################
"""
$(SIGNATURES)
Plot samples of trace on individual subplots.

# Examples
```julia
```

"""
################################################################################
# only implemented for univariate data so far
function plotPredictions(
    trace::Trace{C,A,B},
    transform::TraceTransform,
    data::AbstractArray;
#    model=false,                          # If model <: AbstractModel given, plots true parameter
#    burnin = 0 ,
    CIRegion=[0.025, 0.975],
    dates = false,
#    layout = length_constrained(transform.tagged),
    plotsize=plot_default_size,
#    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {C,A,B}


    # Assign hyperparameter from transform
    burnin = transform.burnin

    #Still contains t+1 predictions at t
    predictions = [trace.diagnostics[chain].base.prediction for chain in eachindex(trace.diagnostics)]
    data_real = data[end-length(predictions)+1:end]
    x_iter_real = dates == false ? collect(1:length(data))[end-length(predictions)+1:end] : dates[end-length(predictions)+1:end]
    ArgCheck.@argcheck length(x_iter_real) == length(data_real) == length(predictions)

    #Remove burnin
    predictions_burnin = predictions[1+burnin:end]
    data_real_burnin = data_real[1+burnin:end]
    x_iter_burnin = x_iter_real[1+burnin:end]

    #Put Data and Prediction on same scale
    predictions_burnin_shifted = predictions_burnin[1:end-1]
    data_real_burnin_shifted = data_real_burnin[2:end]
    x_iter_burnin_shifted = x_iter_burnin[2:end]

    Ndata = length(predictions_burnin_shifted[begin][begin])
    ArgCheck.@argcheck Ndata == 1 "Only univariate prediction plots supported at the moment"

    # Compute Posterior Mean and Credible Interval
    post_mean = zeros(Float64, Ndata, length(predictions_burnin_shifted))
    CI = Matrix{Tuple{Float64,Float64}}(undef, Ndata, length(predictions_burnin_shifted))
    for Nchain in Base.OneTo(length(predictions_burnin_shifted))
        dataᵗ⁺¹ = predictions_burnin_shifted[Nchain]
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

    # Plot data
    plot_data = plot(;
        layout=(Ndata, 1),
        size=plotsize,
        foreground_color_legend = :transparent,
        background_color_legend = :transparent,
        legend=:topleft,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    # Plot CI and data
    for datadim in Base.OneTo(Ndata)
        Plots.plot!(
            x_iter_burnin_shifted, #xiter[end-length(post_mean[datadim, :])+1:end],
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
        Plots.plot!(
            x_iter_burnin_shifted, #xiter[end-Niter+1:end],
            data_real_burnin_shifted, #[start_date:end, datadim];
            linewidth=1.0,
            label="Sample data",
            color="black", #linestyle=:dot,
            subplot=datadim,
        )
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
    # Return Plot
    return plot_data
end


################################################################################
#export
export plotPredictions
