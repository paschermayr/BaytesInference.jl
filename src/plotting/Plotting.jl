"Plotting functions for Baytes"
################################################################################
#Import
include("chains.jl")
include("histograms.jl")
include("predictions.jl")
include("credibleinterval.jl")

#Legacy code
include("legacy/Plotting.jl")

#include("priorprediction.jl")
#include("credibleinterval.jl")

################################################################################
#export
export plot
