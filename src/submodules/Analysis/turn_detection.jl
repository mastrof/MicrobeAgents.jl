export rundurations, mean_runduration, mean_turnrate, detect_turns, has_turned


"""
    rundurations(adf, Δt; threshold_angle=0.0)
Evaluate the duration of all runs observed during a simulation.
The dataframe `adf` should contain the field `:vel`.
`Δt` is the integration timestep of the simulation.
`threshold_angle` defines the threshold (in radians) to define a reorientation.
The function returns a vector of run durations for each microbe.

Also works if `adf` is an `AbstractMatrix`, with unit velocity vectors
for each microbe sorted by column, or if `adf` is an `AbstractVector`
(unit velocity vectors for a single microbe).
"""
function rundurations(adf, Δt; threshold_angle=0.0)
    turns = detect_turns(adf; threshold_angle)
    map(col -> diff([0; findall(col)]) .* Δt, eachcol(turns))
end


"""
    mean_runduration(adf, Δt; threshold_angle=0.0)
Evaluate the mean run duration sampled by the simulation.
The dataframe `adf` should contain the field `:vel`.
`Δt` is the integration timestep of the simulation.
`threshold_angle` defines the threshold (in radians) to define a reorientation.

Also works if `adf` is an `AbstractMatrix`, with unit velocity vectors
for each microbe sorted by column, or if `adf` is an `AbstractVector`
(unit velocity vectors for a single microbe).
"""
mean_runduration(adf, Δt; threshold_angle=0.0) = 1.0 / mean_turnrate(adf, Δt; threshold_angle)

"""
    mean_turnrate(adf, Δt; threshold_angle=0.0)
Evaluate the mean turn rate sampled by the simulation.
The dataframe `adf` should contain the field `:vel`.
`Δt` is the integration timestep of the simulation.
`threshold_angle` defines the threshold (in radians) to define a reorientation.

Also works if `adf` is an `AbstractMatrix`, with unit velocity vectors
for each microbe sorted by column, or if `adf` is an `AbstractVector`
(unit velocity vectors for a single microbe).
"""
function mean_turnrate(adf, Δt; threshold_angle=0.0)
    turns = detect_turns(adf; threshold_angle)
    count(turns) / (length(turns) * Δt)
end

"""
    detect_turns(adf; threshold_angle=0.0)
Detect reorientations in the microbe trajectories in `adf`.
Requires `adf` to have a `:vel` field containing the unit-norm velocity vector.
In order to be detected, a reorientation must be larger than the `threshold_angle`
(in radians); the threshold defaults to 0.

Also works if `adf` is an `AbstractMatrix`, with unit velocity vectors
for each microbe sorted by column, or if `adf` is an `AbstractVector`
(unit velocity vectors for a single microbe).
"""
function detect_turns(adf; threshold_angle=0.0)
    vel = vectorize_adf_measurement(adf, :vel)
    detect_turns(vel; threshold_angle)
end

function detect_turns(vel::AbstractMatrix; threshold_angle=0.0)
    mapslices(u -> detect_turns(u; threshold_angle), vel, dims=1)
end

function detect_turns(vel::AbstractVector; threshold_angle=0.0)
    a₀ = UnitRange(1, length(vel)-1)
    a₁ = UnitRange(2, length(vel))
    has_turned.(view(vel,a₀), view(vel,a₁); threshold_angle)
end

"""
    has_turned(v₁, v₂; threshold_angle=0.0)
Determine if `v₂` is rotated w.r.t `v₁` by comparing the angle between
the two vectors with a `threshold_angle` (defaults to 0).
"""
has_turned(v₁, v₂; threshold_angle=0.0) = safe_acos(dot(v₁,v₂)) > abs(threshold_angle)%π

"""
    safe_acos(x)
Numerically safe version of `acos`.
Useful when `x` exceeds [-1,1] due to floating point errors.
"""
safe_acos(x::T) where T = x ≈ 1 ? zero(x) : x ≈ -1 ? T(π) : acos(x)