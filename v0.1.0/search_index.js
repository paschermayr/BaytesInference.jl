var documenterSearchIndex = {"docs":
[{"location":"intro/#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"intro/","page":"Introduction","title":"Introduction","text":"Yet to be done.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = BaytesInference","category":"page"},{"location":"#BaytesInference","page":"Home","title":"BaytesInference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BaytesInference.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BaytesInference]","category":"page"},{"location":"#BaytesInference.BaytesInference","page":"Home","title":"BaytesInference.BaytesInference","text":"Bayesian inference on state space models\n\n\n\n\n\n","category":"module"},{"location":"#BaytesInference.check_pf-Union{Tuple{T}, Tuple{Random.AbstractRNG, BaytesFilters.ParticleFilter, ModelWrappers.Objective, T}} where T<:NamedTuple","page":"Home","title":"BaytesInference.check_pf","text":"check_pf(_rng, pf, objective, θ; trials, exact, dmax, plotsize, param_color)\n\n\nCheck log target estimate variance for a interval of different model parameter values.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.compute_escore-Union{Tuple{R}, Tuple{Vector{R}, R}} where R<:Real","page":"Home","title":"BaytesInference.compute_escore","text":"compute_escore(prediction, dataₜ)\n\n\nCompute escore of trace given data.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.compute_ℓscore-Union{Tuple{D}, Tuple{T}, Tuple{S}, Tuple{Array{Vector{S}, 1}, Array{Vector{T}, 1}, D}} where {S, T<:Real, D}","page":"Home","title":"BaytesInference.compute_ℓscore","text":"compute_ℓscore(predictionᵛ, ℓweightsₙ, data; burnin)\n\n\nCompute log score of trace given data.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.flatten_chains-Union{Tuple{F}, Tuple{Baytes.Trace, ModelWrappers.ModelWrapper, Integer, F}} where F<:ModelWrappers.FlattenDefault","page":"Home","title":"BaytesInference.flatten_chains","text":"Same as above, but only flatten parameter that have valid constraint.\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.flatten_index-Union{Tuple{F}, Tuple{Baytes.Trace, ModelWrappers.ModelWrapper, Integer, F}} where F<:ModelWrappers.FlattenDefault","page":"Home","title":"BaytesInference.flatten_index","text":"Flatten subset of model.val for each chain at index 'index'. Useful if indices at different chains are compared.\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.get_vals-Tuple{Baytes.Trace, Integer}","page":"Home","title":"BaytesInference.get_vals","text":"Return all parameter values, removing first burnin draws.\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.hpd_grid","page":"Home","title":"BaytesInference.hpd_grid","text":"hpd_grid(samples)\nhpd_grid(samples, alpha; iter, roundto, ϵ)\n\n\nCompute HPD interval of a vector of MCMC samples.\n\nExamples\n\n\n\n\n\n\n\n","category":"function"},{"location":"#BaytesInference.plotChain-Tuple{Any, Any}","page":"Home","title":"BaytesInference.plotChain","text":"plotChain(trace, algorithm)\n\n\nPlot samples of trace.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.plotDiagnostics-Tuple{Any}","page":"Home","title":"BaytesInference.plotDiagnostics","text":"plotDiagnostics(diagnostics)\n\n\nPlot diagnostics of trace.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.plotLatent-Tuple{Any, Any}","page":"Home","title":"BaytesInference.plotLatent","text":"plotLatent(trace, algorithm)\n\n\nPlot latent state trajectory of trace. Only supported for state space models with HMM or HSMM transitions.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.plotPosteriorPrediction-Tuple{Any, Any}","page":"Home","title":"BaytesInference.plotPosteriorPrediction","text":"plotPosteriorPrediction(trace, algorithm)\n\n\nPlot posterior predictive samples of trace. Only supported for some algorithms.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.plot_escore-Union{Tuple{D}, Tuple{Baytes.Trace, Any, D}} where D","page":"Home","title":"BaytesInference.plot_escore","text":"Plot escore of trace given data.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"},{"location":"#BaytesInference.plot_logscore-Union{Tuple{D}, Tuple{Baytes.Trace, Any, D}} where D","page":"Home","title":"BaytesInference.plot_logscore","text":"plot_logscore(trace, smc, data; burnin, delete∞, plotsize, param_color, fontsize, axissize)\n\n\nPlot log score of trace given data.\n\nExamples\n\n\n\n\n\n\n\n","category":"method"}]
}