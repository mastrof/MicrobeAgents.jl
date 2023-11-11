export acf!, acf

function Autocorrelations.acf!(df::AbstractDataFrame, sym; kwargs...)
    gdf = groupby(df, :id)
    acf!(gdf, sym; kwargs...)
end
function Autocorrelations.acf!(df::GroupedDataFrame, sym; kwargs...)
    f(x) = acf(x; kwargs...)
    transform!(df, sym => f => Symbol("$(sym)_acf"))
end
function Autocorrelations.acf(df::AbstractDataFrame, sym; kwargs...)
    acf(groupby(df, :id), sym; kwargs...)
end
function Autocorrelations.acf(df::GroupedDataFrame, sym; kwargs...)
    f(x) = acf(x; kwargs...)
    [f(g[!,sym]) for g in df]
end
