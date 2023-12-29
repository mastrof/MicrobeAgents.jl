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
    nsteps = unique(adf[!,:step]) |> length
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
function MeanSquaredDisplacement.unfold!(adf::AbstractDataFrame, model::ABM; key=:position)
    gdf = groupby(adf, :id)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(adf::AbstractDataFrame, extent; key=:position)
    gdf = groupby(adf, :id)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(gdf::GroupedDataFrame, model::ABM; key=:position)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(gdf::GroupedDataFrame, extent; key=:position)
    newkey = Symbol(join((string(key), "unfold"), '_'))
    transform!(gdf, key => (x -> unfold(x, extent)) => newkey)
end
