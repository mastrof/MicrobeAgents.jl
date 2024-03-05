export driftvelocity_point!, driftvelocity_direction!

"""
    driftvelocity_point!(adf, target; normalize=false, kwargs...)
    driftvelocity_point!(gdf, target; normalize=false, kwargs...)
Evaluate the drift velocity of microbes towards a `target` point,
extracting their positions and velocities from the agent dataframe `adf`
or from the grouped dataframe `gdf` (grouped by `:id`).
Requires the existence of columns containing the position and velocity
of microbes.
By convention the drift velocity will be positive
for motion towards the target point.

**Keywords**
- `normalize = false`: whether drift velocities should be normalized by
  the instantaneous speed of the microbes
- `pos_key::Symbol = :position`: name of the column containing microbe positions
- `vel_key::Symbol = :velocity`: name of the column containing microbe velocities
"""
function driftvelocity_point! end

function driftvelocity_point!(
    adf::AbstractDataFrame,
    target::Union{<:SVector{D},<:NTuple{D}};
    pos_key::Symbol=:position,
    vel_key::Symbol=:velocity,
    normalize=false,
) where D
    gdf = groupby(adf, :id)
    driftvelocity_direction!(gdf, target; pos_key, vel_key, normalize)
end
function driftvelocity_point!(
    gdf::GroupedDataFrame,
    target::Union{<:SVector{D},<:NTuple{D}};
    pos_key::Symbol=:position,
    vel_key::Symbol=:velocity,
    normalize=false,
) where D
    df = parent(gdf)
    @assert hasproperty(df, pos_key)
    @assert hasproperty(df, vel_key)
    cols = [pos_key, vel_key]
    drift(pos, vel) = driftvelocity_point(pos, vel, target; normalize)
    transform!(gdf, cols => drift => :drift_point)
end

function driftvelocity_point(
    pos::AbstractVector{<:SVector{D}},
    vel::AbstractVector{<:SVector{D}},
    target::Union{<:SVector{D},<:NTuple{D}};
    normalize=false
) where D
    @assert length(pos) == length(vel)
    v_drift = zeros(length(pos))
    for k in eachindex(pos)
        # vector from microbe to target
        axis = target .- pos[k]
        # project velocity onto axis
        v_drift[k] = dot(vel[k], axis) / norm(axis)
        if normalize
            v_drift[k] /= norm(vel[k])
        end
    end
    return v_drift
end


"""
    driftvelocity_direction!(adf, target; normalize=false, kwargs...)
    driftvelocity_direction!(gdf, target; normalize=false, kwargs...)
Evaluate the drift velocity of microbes along a `target` direction,
extracting their positions and velocities from the agent dataframe `adf`
or from the grouped dataframe `gdf` (grouped by `:id`).
Requires the existence of a column containing the velocity of microbes.
By convention the drift velocity will be positive
for motion along the target directoin.

**Keywords**
- `normalize = false`: whether drift velocities should be normalized by the
  instantaneous speed of the microbes
- `vel_key::Symbol = :velocity`: name of the column containing microbe velocities
"""
function driftvelocity_direction! end

function driftvelocity_direction!(
    adf::AbstractDataFrame,
    target::Union{<:SVector{D},<:NTuple{D}};
    vel_key::Symbol=:velocity,
    normalize=false
) where D
    gdf = groupby(adf, :id)
    driftvelocity_direction!(gdf, target; vel_key, normalize)
end
function driftvelocity_direction!(
    gdf::GroupedDataFrame,
    target::Union{<:SVector{D},<:NTuple{D}};
    vel_key::Symbol=:velocity,
    normalize=false
) where D
    df = parent(gdf)
    @assert hasproperty(df, vel_key)
    drift(vel) = driftvelocity_direction(vel, target; normalize)
    transform!(gdf, vel_key => drift => :drift_direction)
end

function driftvelocity_direction(
    vel::AbstractVector{<:SVector{D}},
    target::Union{<:SVector{D},<:NTuple{D}};
    normalize=false
) where D
    if normalize
        v_drift = map(v -> dot(v, target) / (norm(target)*norm(v)), vel)
    else
        v_drift = map(v -> dot(v, target) / norm(target), vel)
    end
    return v_drift
end
