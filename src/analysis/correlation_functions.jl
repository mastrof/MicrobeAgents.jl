export autocorrelation

"""
    autocorrelation(adf, sym::Union{String,Symbol}; tstep::Int=1)
Return the non-normalized autocorrelation function of quantity `sym` extracted
from the dataframe `adf` (`adf[!,sym]`).
The kwarg `tstep` defines the values of time-lags for sampling (defaults to 1).
"""
function autocorrelation(adf, sym::Union{String,Symbol}; tstep::Int=1)
    timeseries = vectorize_adf_measurement(adf, sym)
    autocorrelation(timeseries; tstep=tstep)
end # function

function autocorrelation(timeseries::T; tstep::Int=1) where {S,T<:AbstractVector{S}}
    nsteps, = size(timeseries)
    timelags = range(0, nsteps-2; step=tstep)
    f = zeros(nsteps-1)
    for (j,Δt) in enumerate(timelags)
        for t₀ in 1:nsteps-Δt
            u = timeseries[t₀]
            v = timeseries[t₀+Δt]
            f[j] += dot(u,v)
        end # for
        f[j] /= (nsteps-j)
    end # for
    return f
end # function

function autocorrelation(timeseries::T; tstep::Int=1) where {S,T<:AbstractMatrix{S}}
    nmicrobes, nsteps = size(timeseries)
    timelags = range(0, nsteps-2; step=tstep)
    f = zeros(nsteps-1)
    for (j,Δt) in enumerate(timelags)
        for t₀ in 1:nsteps-Δt
            for i in 1:nmicrobes
                u = timeseries[i, t₀]
                v = timeseries[i, t₀+Δt]
                f[j] += dot(u,v)
            end # for
        end # for
        f[j] /= (nmicrobes * (nsteps-j))
    end # for
    return f
end # function