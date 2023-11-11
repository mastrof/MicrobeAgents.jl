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

function MeanSquaredDisplacement.unfold!(adf::AbstractDataFrame, model::ABM)
    gdf = groupby(adf, :id)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(adf::AbstractDataFrame, extent::SVector{D}) where {D}
    gdf = groupby(adf, :id)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(gdf::GroupedDataFrame, model::ABM)
    extent = spacesize(model)
    unfold!(gdf, extent)
end
function MeanSquaredDisplacement.unfold!(gdf::GroupedDataFrame, extent::SVector{D}) where {D}
    transform!(gdf, :pos => (x -> unfold(x, extent)) => :pos_unfold)
end
