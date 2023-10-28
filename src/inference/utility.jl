
################################################################################
# Get MCMC Chain from Trace and Transform

"""
$(SIGNATURES)
Change trace to ArviZ Inference Object.

# Examples
```julia
```

"""
function to_ArviZ end

"""
$(SIGNATURES)
Change trace to MCMCChain.

# Examples
```julia
```

"""
function to_MCMCChain end

function trace_to_MCMCChain(trace::Trace, transform::TraceTransform)
    arr3D = Baytes.trace_to_3DArray(trace, transform)
#    mcmcchain = MCMCChains.Chains( permutedims(arr3D, (1, 3, 2) ), transform.paramnames )
    mcmcchain = to_MCMCChain(permutedims(arr3D, (1, 3, 2) ), transform.paramnames)
    return mcmcchain
end
function trace_to_ArviZ(trace::Trace, transform::TraceTransform)
    mcmcchain = trace_to_MCMCChain(trace, transform)
    return to_ArviZ(mcmcchain)
end

################################################################################
#export
export trace_to_MCMCChain, trace_to_ArviZ