export Microbe

"""
    Microbe{D} <: AbstractMicrobe{D}
Basic microbe type for simple simulations.

Default parameters:
- `id::Int = rand(1:typemax(Int))` identifier used internally by Agents.jl
- `pos::NTuple{D,Float64} = ntuple(zero,D)` position
- `motility = RunTumble()` motile pattern
- `vel::NTuple{D,Float64} = rand_vel(D, motility)` velocity vector
- `turn_rate::Float64 = 1.0` frequency of reorientations
- `rotational_diffusivity::Real` coefficient of brownian rotational diffusion
- `radius::Float64 = 0.0` equivalent spherical radius of the microbe
- `state::Float64 = 0.0` generic variable for a scalar internal state
"""
mutable struct Microbe{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility 
    vel::NTuple{D,Float64}
    turn_rate::Float64
    rotational_diffusivity::Float64
    radius::Float64
    state::Float64

    Microbe{D}(
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility::AbstractMotility = RunTumble(),
        vel::NTuple{D,<:Real} = rand_vel(D, motility),
        turn_rate::Real = 1.0,
        rotational_diffusivity::Real = 0.0,
        radius::Real = 0.0,
        state::Real = 0.0,
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(turn_rate),
        Float64(rotational_diffusivity), Float64(radius), Float64(state)
    )
end # struct

r2dig(x) = round(x, digits=2)
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D}) where D
    println(io, "$(typeof(m)) with $(typeof(m.motility)) motility")
    println(io, "position (μm): $(r2dig.(m.pos)); velocity (μm/s): $(r2dig.(m.vel))")
    println(io, "average unbiased turn rate (Hz): $(r2dig(m.turn_rate))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate])
    print(io, "other properties: " * join(s, ", "))
end