################################################################################
"""
$(SIGNATURES)
Compute WAIC. 'incremental_ℓlikᵥ' is a vector of vector, where
the inner vector contains the incremental likelihood at each data point,
and the outer vector has length equal to the number of MCMC iterations.

# Examples
```julia
```

"""
function compute_waic(incremental_ℓlikᵥ::Vector{Vector{F}}, t₀::Integer = 1) where {F<:Real}
        #!NOTE: Inner vector is of length T, where T == # of data points.
        #!NOTE2: Outter vector is # of MCMC samples
    T = length(incremental_ℓlikᵥ[begin])
    N = length(incremental_ℓlikᵥ)
    # Compute lppd
    lppdᵥ = [ BaytesCore.logmeanexp( [ incremental_ℓlikᵥ[iter][t] for iter in Base.OneTo(N)] ) for t in t₀:T ]
    lppdᵥ_nan = findall(isnan.(lppdᵥ))
    if length(lppdᵥ_nan) > 0
        println("NANs found for lppdᵥ: ", length(lppdᵥ_nan), " at outter index ", lppdᵥ_nan)
    end
    lppd = sum( filter(!isnan, lppdᵥ) )
    # Compute p_waic
    p_waicᵥ = [ var(incremental_ℓlikᵥ[iter][t] for iter in Base.OneTo(N); corrected=false) for t in t₀:T ]
    p_waicᵥ_nan = findall(isnan.(p_waicᵥ))
    if length(p_waicᵥ_nan) > 0
        println("NANs found for p_waicᵥ: ", length(p_waicᵥ_nan), " at outter index ", p_waicᵥ_nan)
    end
    p_waic = sum( filter(!isnan, p_waicᵥ) )
    # Compute WAIC
    waic = -2*lppd +2*p_waic
    # Return output - The LOWER waic, the better. p_waic is a correction term for additional parameter (the higher the more parameter)
    (waic = waic, lppd = lppd, p_waic = p_waic)
end

################################################################################
#export
    export compute_waic
