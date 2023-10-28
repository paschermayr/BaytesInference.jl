module BaytesInferenceArviZExt

############################################################################################
#using BaytesOptim
import BaytesInference: BaytesInference, to_ArviZ, to_MCMCChain
using MCMCChains, ArviZ


function to_MCMCChain(array, names)
    return MCMCChains.Chains( array, names )
end
function to_ArviZ(mcmcchain)
    ArviZ_inference_data = ArviZ.from_mcmcchains(
        mcmcchain
    )
    return ArviZ_inference_data
end

#=
# Get MCMC Chain from Trace and Transform
function trace_to_MCMCChain(trace::Trace, transform::TraceTransform)
    arr3D = Baytes.trace_to_3DArray(trace, transform)
    mcmcchain = MCMCChains.Chains( permutedims(arr3D, (1, 3, 2) ), transform.paramnames )
    return mcmcchain
end
function trace_to_InferenceData(trace::Trace, transform::TraceTransform)
    mcmcchain = trace_to_MCMCChain(trace, transform)
    ArviZ_inference_data = from_mcmcchains(
        mcmcchain
    )
    return ArviZ_inference_data
end
=#

end