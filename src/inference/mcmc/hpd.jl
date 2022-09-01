################################################################################
#
function get_intervals(x::AbstractVector)
    @argcheck length(x) >= 1
    #even index values
    x_even = view(x, 2:2:length(x))
    #odd index values
    x_odd = view(x, 1:2:length(x))
    return [(x_odd[iter], x_even[iter]) for iter in eachindex(x_even)]
end

#=
NOTE: Function taken and adjusted from Python file:
https://github.com/aloctavodia/BAP/blob/master/first_edition/code/Chp1/hpd.py
=#
"""
$(SIGNATURES)
Compute HPD interval of a vector of MCMC samples.

# Examples
```julia
```

"""
function hpd_grid(samples, alpha=0.05;
    iter = 10^4,
    roundto = 2,
    ϵ = 1.0e-3
)
    # get upper and lower bounds
    l = minimum(samples)
    u = maximum(samples)
    x = LinRange(l, u, iter)
    # Get KDE estimate
    #density = KernelDensity.kde(samples; kernel=Normal)
    density = KernelDensity.kde(samples)
    y = KernelDensity.pdf(density, x)
    # Compute normalized density, and return tuple of it along x range
    sum_y = sum(y)
    y = y ./ sum_y
    # Sort x and normalizd y based on normalized y from highest to lowest
    idx_sorted = sortperm(y, rev=true)
    x = x[idx_sorted]
    y = y[idx_sorted]
    #Set hdv and cum sum and loop through until quantile reach
    y_cum_sum = 0.0
    hdv = Float64[]
    for (index, val) in enumerate(y)
        y_cum_sum += val
        push!(hdv, x[index])
        if y_cum_sum >= (1-alpha)
            break
        end
    end
    #Sort quantiles
    hdv = sort(hdv)
    diff = (u-l)/50  # differences of 5%
    hpd = Float64[]
    push!(hpd, round(minimum(hdv); digits = roundto))
    #
    for index in 2:length(hdv)
        if hdv[index]-hdv[index-1] >= diff
            push!(hpd, round(hdv[index-1]; digits = roundto))
            push!(hpd, round(hdv[index]; digits = roundto))
        end
    end
    push!(hpd, round(maximum(hdv); digits = roundto))
    hpd = get_intervals(hpd)
    modes = Float64[]
    for value in hpd
         x_hpd = x[(x .> value[1] - ϵ) .& (x .< value[2] + ϵ)]
         y_hpd = y[(x .> value[1] - ϵ) .& (x .< value[2] + ϵ)]
         push!(modes, round(x_hpd[argmax(y_hpd)]; digits = roundto))
    end
    return hpd, x, y, modes
end

################################################################################
#export
    export hpd_grid


#=
using KernelDensity, Distributions, ArgCheck, Plots

d = MixtureModel([Normal(-10, 2), Normal(0, 2), Normal(10, 2)], [.33, .33, .34])
iter = 10^4
samples = rand(d, iter)

_alpha = 0.1
hpd, x, y, modes = hpd_grid(samples, _alpha)

_fonttemp = 12
p = plot(;
    layout=(1, 1),
    size=(1000,1000),
    legend=:topleft,
    xguidefontsize=_fonttemp,
    yguidefontsize=_fonttemp,
    legendfontsize=_fonttemp,
    xtickfontsize=_fonttemp,
    ytickfontsize=_fonttemp,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
)
histogram!(samples, bins=100, alpha = 0.5)
for iter in eachindex(hpd)
    vspan!([hpd[iter][1], hpd[iter][2]], linecolor = :grey, fillcolor = :grey, alpha = 0.25, label=string("Mode ", iter))
end
#Check against quantiles
quant1 = quantile(samples, _alpha/2)
quant2 = quantile(samples, 1 - _alpha/2)
vline!([quant1], color="black", label="2.5%CI")
vline!([quant2], color="black", label="97.5%CI")
xaxis!("90%CI and HPD for a Normal Mixture model with 3 components.")
p
=#
