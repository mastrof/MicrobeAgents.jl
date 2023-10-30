export acf!

function Autocorrelations.acf!(df::GroupedDataFrame, sym; kwargs...)
    f(x) = acf(x; kwargs...)
    transform!(df, sym => f => Symbol("$(sym)_acf"))
end
