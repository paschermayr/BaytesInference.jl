################################################################################
"""
$(SIGNATURES)
Compute DIC. 'ℓlikθ' is the log target evaluation given posterior mean.'ℓlikᵥ' is a vector of target evaluations p(data | θ_n) for n=1:N. Note that target in this case is likelihood, not posterior.

# Examples
```julia
```

"""
function compute_dic(ℓlikθ::F, ℓlikᵥ::Vector{F}, t₀::Integer = 1) where {F<:Real}
    ℓlikᵥ_nan = findall(isnan.(ℓlikᵥ))
    if length(ℓlikᵥ_nan) > 0
        println("NANs found for ℓlikᵥ: ", length(ℓlikᵥ_nan), " at outter index ", ℓlikᵥ_nan)
    end
    p_dic = 2 * ( ℓlikθ - mean(ℓlikᵥ) )
    dic = - 2 * ℓlikθ + 2 * p_dic
    # Return output - The LOWER dic, the better. p_DIC is a correction term for additional parameter (the higher the more parameter)
    (dic = dic, ℓlikθ = ℓlikθ, p_dic = p_dic)
end

################################################################################
#export
    export compute_dic
