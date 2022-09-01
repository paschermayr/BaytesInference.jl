"Bayesian inference on state space models"
module BaytesInference

################################################################################
#Import modules
using ModelWrappers
using BaytesCore, BaytesMCMC, BaytesFilters, BaytesPMCMC, BaytesSMC, Baytes

#Utility tools
using DocStringExtensions:
    DocStringExtensions, TYPEDEF, TYPEDFIELDS, FIELDS, SIGNATURES, FUNCTIONNAME
using ArgCheck: ArgCheck, @argcheck
using UnPack: UnPack, @unpack, @pack!
using Random: Random, AbstractRNG, GLOBAL_RNG, randexp
using Statistics: Statistics, mean, std, sqrt, quantile, var
using ProgressMeter

using Plots
import Plots: Plots, plot

#Inference modules
using Distributions
using KernelDensity: KernelDensity, kde, pdf

################################################################################
# Define constants
const plot_default_color = :rainbow_bgyrm_35_85_c71_n256
const plot_default_size = (1000, 1000)
const _fontsize = 12
const _axissize = 12

################################################################################
#Helper structs so we can make plotting for SMC easier
struct _IBIS end
struct _SMC2 end

################################################################################
#Abstract types to be dispatched in Examples section
include("plotting/Plotting.jl")
include("inference/Inference.jl")

################################################################################
export
#Baytes
    plot,
    _SMC2,
    _IBIS

end