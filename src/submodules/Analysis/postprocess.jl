export adf_to_matrix, adf_to_vectors, unfold, unfold!

"""
    adf_to_matrix(adf, sym)
Collect quantity `sym` from the agent dataframe `adf` and return it
in matrix form with dimensions (times, microbes).
Requires all microbes to exist from start to end of the simulations 
(no removals or insertions).
"""
function adf_to_matrix(adf, sym)
    nmicrobes = unique(adf[!,:id]) |> length
    nsteps = unique(adf[!,:time]) |> length
    datatype = typeof(adf[1,sym])
    s = Matrix{datatype}(undef, nsteps, nmicrobes)
    for i in 1:nmicrobes
        for t in 1:nsteps
            s[t,i] = adf[i + (t-1)*nmicrobes, sym]
        end
    end
    return s
end

"""
    adf_to_vectors(adf, sym)
Collect quantity `sym` from the agent dataframe `adf` and return a
vector of such quantity for each microbe id.
"""
function adf_to_vectors(adf, sym)
    gdf = groupby(adf, :id)
    s = [g[!,sym] for g in gdf]
    return s
end


"""
    unfold!(adf, model; key=:position)
    unfold!(adf, extent; key=:position)
    unfold!(gdf, model; key=:position)
    unfold!(gdf, extent; key=:position)
Unfold the periodic trajectories contained in the agent dataframe `adf` or in the grouped
dataframe `gdf` where each microbe trajectory was grouped by `:id`.
The second argument can be either the spatial extent of the model, or the model itself
from which the extent is automatically extracted.

The keyword argument `key=:position` determines what column of the dataframe has to be
unfolded.
"""
function unfold!(adf::AbstractDataFrame, model::ABM; key=:position)
    gdf = groupby(adf, :id)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function unfold!(adf::AbstractDataFrame, extent; key=:position)
    gdf = groupby(adf, :id)
    unfold!(gdf, extent)
end
function unfold!(gdf::GroupedDataFrame, model::ABM; key=:position)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function unfold!(gdf::GroupedDataFrame, extent; key=:position)
    newkey = Symbol(join((string(key), "unfold"), '_'))
    transform!(gdf, key => (x -> unfold(x, extent)) => newkey)
end

"""
    unfold(x::AbstractVector, P)
Unfold timeseries `x` from a periodic domain extending from `0` to `P`.
"""
function unfold(x::AbstractVector, P)
    y = deepcopy(x)
    unfold!(y, P)
    return y
end

"""
    unfold!(x::AbstractVector, P)
Unfold timeseries `x` from a periodic domain extending from `0` to `P`.
The periodicity `P` can be either a real number (for a cubic domain) or
a collection (`AbstractVector` or `NTuple`) with one value for each dimension.
"""
function unfold!(x::AbstractVector, P)
    indices = eachindex(x)
    ind_prev = @view indices[1:end-1]
    ind_next = @view indices[2:end]
    for (i,j) in zip(ind_prev, ind_next)
        x0 = x[i]
        x1 = x[j]
        x[j] = unfold(x1, x0, P)
    end
    return x
end

NTupleOrVec = Union{NTuple{D,<:Real},AbstractVector{<:Real}} where D
function unfold(x1::AbstractVector, x0::AbstractVector, P::Real)
    @assert length(x1) == length(x0)
    map(i -> unfold(x1[i], x0[i], P), eachindex(x1))
end
function unfold(x1::NTuple{D}, x0::NTuple{D}, P::Real) where D
    ntuple(i -> unfold(x1[i], x0[i], P), D)
end
function unfold(x1::AbstractVector, x0::AbstractVector, P::NTupleOrVec)
    @assert length(x1) == length(x0) == length(P)
    map(i -> unfold(x1[i], x0[i], P[i]), eachindex(x1))
end
function unfold(x1::NTuple{D}, x0::NTuple{D}, P::NTupleOrVec) where D
    @assert D == length(P)
    ntuple(i -> unfold(x1[i], x0[i], P[i]), D)
end

"""
    unfold(x1::Real, x0::Real, P::Real)
Unfold the value of `x1` with respect to `x0` from a
domain of periodicity `P`.
"""
function unfold(x1::Real, x0::Real, P::Real)
    dx = x1 - x0
    a = round(abs(dx / P))
    if abs(dx) > P/2
        return x1 - a*P*sign(dx)
    else
        return x1
    end
end

