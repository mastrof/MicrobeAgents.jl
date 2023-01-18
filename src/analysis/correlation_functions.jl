export acf

"""
    acf(x, sym, lags=0:last(x.step)-1)
Evaluate the autocorrelation function for quantity `sym` from dataframe `x` (`x.sym`)
at timelags `lags`.
If `lags` is unspecified, the acf is evaluated for all available timepoints.
"""
function acf(x, sym, lags=0:last(x.step))
    y = vectorize_adf_measurement(x, sym)
    acf(y, lags)
end

"""
    acf(x, lags=0:size(x,1)-1)
Evaluate the autocorrelation function (acf) of a signal `x` at timelags `lags`.
If `lags` is unspecified, the acf is evaluated for all available time points.

The signal `x` can be an `AbstractVector` or an `AbstractMatrix`.
In the latter case, it is assumed that each column represents a distinct signal
and each row a timepoint. The acf is then evaluated for each individual
timeseries and the results are averaged together.

The result is *not* normalized. If `f = acf(x)`, just evaluate
`f./f[1]` to normalize it.
"""
function acf(x::AbstractMatrix, lags=0:size(x,1)-1)
    mean(mapslices(y -> acf(y, lags), x, dims=1), dims=2)
end
function acf(x::AbstractVector{<:Tuple}, lags=0:size(x,1)-1)
    l = size(x,1)
    y = [getindex.(x,i) for i in eachindex(first(x))]
    A = sum([conv(s, reverse(s)) for s in y])
    [A[k+l]/(l-k) for k in lags]
end
function acf(x::AbstractVector{<:Number},lags=0:size(x,1)-1)
    l = size(x,1)
    A = conv(x, reverse(x))
    [A[k+l]/(l-k) for k in lags]
end

