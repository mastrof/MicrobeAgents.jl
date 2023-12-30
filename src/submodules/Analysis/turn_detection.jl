export detect_turns!, run_durations


"""
    run_durations(adf; turns_key=:has_turned)
Evaluate the duration (in frames) of all runs detected during a simulation.
The dataframe should already contain a column of Bools with values `true` for
frames when a turn has occurred and `false` otherwise (see `detect_turns!`).
The function returns a vector of run durations for each microbe.

`turns_key` defines the name of the column where the turn statistics is
stored (defaults to `:has_turned`).
"""
function run_durations(adf::AbstractDataFrame; turns_key::Symbol=:has_turned)
    run_durations(groupby(adf, :id); turns_key)
end
function run_durations(gdf::GroupedDataFrame; turns_key::Symbol=:has_turned)
    [diff([0; findall(getproperty(g, turns_key))]) for g in gdf]
end


"""
    detect_turns!(adf; threshold_angle=0.0, kwargs...)
Detect reorientations in the microbe trajectories in `adf`.
Requires `adf` to have a velocity column.
In order to be detected, a reorientation must be larger than the `threshold_angle`
(in radians); the threshold defaults to 0.

A new column is added to the dataframe, with values `true` when the microbe
has turned with respect to the previous frame, or `false` otherwise.

**Keywords**
- `vel_key::Symbol = :velocity`: name of the column containing microbe velocities
- `new_key::Symbol = :has_turned`: name of the new column storing turn statistics
"""
function detect_turns! end

function detect_turns!(adf::AbstractDataFrame; threshold_angle=0.0,
    vel_key::Symbol=:velocity,
    new_key::Symbol=:has_turned,
)
    detect_turns!(groupby(adf, :id); threshold_angle, vel_key, new_key)
end
function detect_turns!(gdf::GroupedDataFrame; threshold_angle=0.0,
    vel_key::Symbol=:velocity,
    new_key::Symbol=:has_turned,
)
    @assert hasproperty(parent(gdf), vel_key)
    detect(vel) = detect_turns(vel; threshold_angle)
    transform!(gdf, vel_key => detect => new_key)
end
function detect_turns(vel::AbstractVector{<:SVector}; threshold_angle=0.0)
    l = length(vel)
    i0 = UnitRange(1, l-1)
    i1 = UnitRange(2, l)
    [false; has_turned.(view(vel, i0), view(vel, i1); threshold_angle)]
end

"""
    has_turned(v₁, v₂; threshold_angle=0.0)
Determine if `v₂` is rotated w.r.t `v₁` by comparing the angle between
the two vectors with a `threshold_angle` (defaults to 0).
"""
function has_turned(v₁, v₂; threshold_angle=0.0)
    u = dot(v₁,v₂) / (norm(v₁)*norm(v₂))
    safe_acos(u) > abs(threshold_angle)%π
end

"""
    safe_acos(x)
Numerically safe version of `acos`.
Useful when `x` exceeds [-1,1] due to floating point errors.
"""
safe_acos(x::T) where T = x ≈ 1 ? zero(x) : x ≈ -1 ? T(π) : acos(x)
