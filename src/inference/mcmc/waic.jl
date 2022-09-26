################################################################################
#
"""
$(SIGNATURES)
Compute WAIC. 'incremental_ℓlikᵥ' is a vector of vector, where the inner vector contains the incremental likelihood at each data point, and the outer vector has length equal to the number of MCMC iterations.

# Examples
```julia
```

"""
function compute_waic(incremental_ℓlikᵥ::Vector{Vector{F}}, initial::Integer) where {F<:Real}
    T = length(incremental_ℓlikᵥ[begin])
    N = length(incremental_ℓlikᵥ)
    # Compute lppd
    lppd = sum( BaytesCore.logmeanexp( [ incremental_ℓlikᵥ[iter][t] for iter in Base.OneTo(N)] ) for t in t₀:T )
    # Compute p_waic
    p_waic = sum( var(incremental_ℓlikᵥ[iter][t] for iter in Base.OneTo(N); corrected=false) for t in t₀:T )
    # Compute WAIC
    waic = -2*lppd +2*p_waic
    # Return output
    (waic = waic, lppd = lppd, p_waic = p_waic)
end

################################################################################
#export
    export compute_waic
