################################################################################
"""
$(SIGNATURES)
Check log target estimate variance for a interval of different model parameter values.

# Examples
```julia
```

"""
function filter_forward(objective::Objective)
    return error("No analytical formula to compute log target density provided, dispatch filter_forward(myobjective) with return values a (matrix of state probabilities, logtargetdensity).")
end

################################################################################
#import Main.Baytes:check_pf

#=objective = objectiveᵖᶠ
a = .09:.09:.9
θ = (;
    μ = -2.:.41:2.,
    σ = .2:.2:2.,
    #p = [[a[iter], 1-a[iter] ] for iter in eachindex(a)],
)
len = [length(θ[iter]) for iter in eachindex(θ)]
=#
"""
$(SIGNATURES)
Check log target estimate variance for a interval of different model parameter values.

# Examples
```julia
```

"""
function check_pf(
    _rng::Random.AbstractRNG,
    pf::ParticleFilter,
    objective::Objective,
    #A tuple with keys = Modelparameter that contains a tuple of ranges for the first parameter in the field
    θ::T;
    trials=15,
    exact=true,
    dmax=1 * size(objective.data, 1),
    plotsize=plot_default_size,
    param_color=plot_default_color,
) where {T<:NamedTuple}
    ##
    len = [length(θ[iter]) for iter in eachindex(θ)]
    ## Initiate container that holds likelihood estimates
    loglik_pf = [zeros(trials, length(θ[iter])) for iter in eachindex(θ)]
    loglik_exact = [zeros(length(θ[iter])) for iter in eachindex(θ)]
    ## Calculate Likelihood estimate vs exact solution
    for (idx, sym) in enumerate(keys(θ))
        θ_temp = deepcopy(getfield(objective.model.val, sym))
        model_temp = deepcopy(objective.model)
        println("Computing for parameter: ", sym)
        ProgressMeter.@showprogress "Computing..." for iter in eachindex(θ[sym])
            #Set parameter to current range
            θ_temp[1] = getfield(θ, sym)[iter]
            #Set temporary theta to objective
            ModelWrappers.fill!(
                model_temp, BaytesCore.to_NamedTuple((sym,), copy(θ_temp))
            )
            Base.Threads.@threads for trial in 1:trials
                _, diagnostics = propose!(_rng,
                    deepcopy(pf),
                    deepcopy(model_temp),
                    objective.data
                )
                loglik_pf[idx][trial, iter] = diagnostics.base.ℓobjective
            end
            if exact
                objective_temp = Objective(model_temp, objective.data, objective.tagged)
                _, loglik_exact[idx][iter] = filter_forward(objective_temp)
            end
        end
    end
    ## Plot results
    plot_ll = plot(; layout=(size(loglik_pf, 1), 1), size=plot_default_size, legend=false)
    #    counter = 0

    for (idx, sym) in enumerate(keys(θ))
        #        counter += 1
        #Plots.plot!(new, ll_e, label="exact", color="black")
        Plots.scatter!(
            θ[idx],
            loglik_pf[idx]';
            ylabel="PF ℓlik estimate vs exact",
            markerstrokewidth=0.5,
            alpha=1.0,
            shape=:o,
            markerstrokecolor=Plots.palette(param_color, (length(keys(θ)) + 1))[idx],#"grey",
            markersize=3,
            #label="Posterior Mean, 95%CI",
            color=Plots.palette(param_color, (length(keys(θ)) + 1))[idx],
            subplot=idx,
        )
        if exact
            Plots.plot!(
                θ[idx], loglik_exact[idx]; label="exact", color="black", subplot=idx
            )
        end
        xaxis!(string("Parameter range for ", sym); subplot=idx)
    end
    ## Return plot and estimates
    return plot_ll, loglik_pf, loglik_exact
end

################################################################################
#export
export check_pf, filter_forward
