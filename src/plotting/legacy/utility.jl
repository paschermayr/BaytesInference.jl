############################################################################################
# Utility functions to deconstruct traces to arrays for plotting. Will be depreciated for newer Baytes.jl utility functions.

function get_chaindiagnostics(
    trace::Trace, chain::Integer, Nalgorithm::Integer, burnin::Integer
)
    return @view(trace.diagnostics[chain][Nalgorithm][(1 + burnin):end])
end

function get_chainvals(trace::Trace, chain::Integer, burnin::Integer)
    return @view(trace.val[chain][(1 + burnin):end])
end

"Return all parameter values, removing first burnin draws."
function get_vals(trace::Trace, burnin::Integer)
    Nchains = length(trace.val)
    return reduce(
        vcat, @view(trace.val[iter][(1 + burnin):end]) for iter in Base.OneTo(Nchains)
    )
end

function get_vals(trace::Trace, tagged::Tagged, burnin::Integer)
    Nchains = length(trace.val)
    return map(Base.OneTo(Nchains)) do chain
        samples = get_chainvals(trace, chain, burnin)
        map(iter -> BaytesCore.subset(samples[iter], tagged.parameter), eachindex(samples))
    end
end

"Same as above, but only flatten parameter that have valid constraint."
function flatten_chains(
    trace::Trace, model::ModelWrapper, burnin::Integer, flattenparam::F
) where {F<:FlattenDefault}
    return map(Base.OneTo(length(trace.val))) do chain
        samples = get_chainvals(trace, chain, burnin)
        reduce(
            hcat,
            ModelWrappers.flatten(model.info.reconstruct, samples[iter]) for iter in eachindex(samples)
        )
    end
end
"Flatten subset of model.val for each chain at index 'index'. Useful if indices at different chains are compared."
function flatten_index(
    trace::Trace, model::ModelWrapper, index::Integer, paramflatten::F
) where {F<:FlattenDefault}
    Nchains = length(trace.val)
    return reduce(
        hcat,
        ModelWrappers.flatten(model.info.reconstruct, trace.val[chain][index]) for chain in Base.OneTo(Nchains)
    )
end
function flatten_index(
    trace::Trace, tagged::Tagged, index::Integer, paramflatten::F
) where {F<:FlattenDefault}
    Nchains = length(trace.val)
    return reduce(
        hcat,
        flatten(
                tagged.info.reconstruct,
                BaytesCore.subset(trace.val[chain][index], tagged.parameter)
        )  for chain in Base.OneTo(Nchains)
    )
end

################################################################################
#export
#export flatten_index, get_vals 
