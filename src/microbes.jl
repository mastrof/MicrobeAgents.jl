export Microbe

"""
    Microbe{D} <: AbstractMicrobe{D}
Base microbe type for simple simulations.

Default parameters:
- `motility = RunTumble()` motile pattern
- `turn_rate::Float64 = 1.0` frequency of reorientations
- `rotational_diffusivity::Real` coefficient of brownian rotational diffusion
- `radius::Float64 = 0.0` equivalent spherical radius of the microbe
- `state::Float64 = 0.0` generic variable for a scalar internal state

`Microbe` has the additional required fields
- `id::Int` an identifier used internally
- `pos::SVector{D,Float64}` spatial position
- `vel::SVector{D,Float64}` unit velocity vector
- `speed::Float64` magnitude of the velocity vector
"""
@agent struct Microbe{D}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D}
    speed::Float64
    motility = RunTumble()
    turn_rate::Float64 = 1.0
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    state::Float64 = 0.0
end

r2dig(x) = round(x, digits=2)
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D}) where D
    println(io, "$(typeof(m)) with $(typeof(m.motility)) motility")
    println(io, "position (μm): $(r2dig.(m.pos)); velocity (μm/s): $(r2dig.(m.vel.*m.speed))")
    println(io, "average unbiased turn rate (Hz): $(r2dig(m.turn_rate))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate])
    print(io, "other properties: " * join(s, ", "))
end
