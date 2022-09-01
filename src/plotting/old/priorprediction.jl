################################################################################
"""
$(SIGNATURES)
Plot Model prior distributions.

# Examples
```julia
```

"""
function plotPrior(
    model::Model,
    sym::Symbol;
    Maxparam=10,
    Nsamples=10000,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
)
    ## Check if symbol contains only univariate parameter
    Nparam = model.info[sym].size
    ArgCheck.@argcheck typeof(model.info[sym].prior) <:
                       Union{UnivariateDistribution,Vector{<:UnivariateDistribution}} "plotPrior only works for univariate priors"
    ArgCheck.@argcheck Maxparam > Nparam string(
        "Keyword Maxparam set to ",
        Maxparam,
        ", parameter to plot size is ",
        Nparam,
        ", increase Maxparam if you want to proceed.",
    )
    ## Plot defaults
    plot_prior = plot(;
        layout=(Nparam, 1),
        size=plot_default_size,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Plot parameter for each chain
    if Nparam > 1
        for iter in Base.OneTo(Nparam)
            #Sample from prior
            θ = rand(model.info[sym].prior[iter], Nsamples)
            #Plot histogram
            Plots.histogram!(
                θ;
                ylabel=string(sym, " ", iter),
                color=Plots.palette(param_color, (Nparam + 1))[iter],
                subplot=iter,
            )
        end
        xaxis!("Value"; subplot=Nparam)
    else
        #Plot histogram
        Plots.histogram!(
            rand(model.info[sym].prior, Nsamples);
            ylabel=string(sym),
            color=Plots.palette(param_color, (Nparam + 1)),
            subplot=1,
        )
    end
    return plot_prior
end

################################################################################
"""
$(SIGNATURES)
Plot samples of prior predictive distribution.

# Examples
```julia
```

"""
function plotPriorPrediction(
    model::Model,
    data::D,
    sym=keys(model.info);
    Nsamples=10000,
    plotsize=plot_default_size,
    param_color=plot_default_color,
    fontsize=_fontsize,
    axissize=_axissize,
) where {D}
    ## Check if predict method exists
    ArgCheck.@argcheck hasmethod(predict, Tuple{typeof(model.value),typeof(data)}) "No predict function given to your Model - dispatch BaytesMCMC.predict(model::AbstractModel) to your Model"
    ## Preallocate output
    Ndata = size(predict(model.value, data), 1)
    predictions = Vector{typeof(predict(model.value, data))}(undef, Nsamples)
    modelᵗᵉᵐᵖ = deepcopy(model)
    ## Plot defaults
    plot_prior = plot(;
        layout=(Ndata, 1),
        size=plot_default_size,
        legend=false,
        xguidefontsize=fontsize,
        yguidefontsize=fontsize,
        legendfontsize=fontsize,
        xtickfontsize=axissize,
        ytickfontsize=axissize,
    )
    ## Sample new model parameter from prior and predict new data
    for iter in Base.OneTo(Nsamples)
        #Sample new paramter from prior
        BaytesTools.sample!(modelᵗᵉᵐᵖ, sym)
        #Predict new parameter
        predictions[iter] = predict(modelᵗᵉᵐᵖ.value, data)
    end
    ## Plot prior predictive distribution
    for datadim in Base.OneTo(Ndata)
        Plots.histogram!(
            getindex.(predictions, datadim);
            ylabel="Prior predictive distribution",
            color=Plots.palette(param_color, (Ndata + 1))[datadim],
            subplot=datadim,
        )
    end
    xaxis!("Value"; subplot=Ndata)
    ## Return plot
    return plot_prior
end

################################################################################
#export
export plotPrior, plotPriorPrediction
