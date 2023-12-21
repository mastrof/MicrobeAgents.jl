export Microbe

"""
    Microbe{D,M} <: AbstractMicrobe{D,M}
Base microbe type for simple simulations.

`Microbe` has the required fields
- `id::Int` an identifier used internally
- `pos::SVector{D,Float64}` spatial position
- `vel::SVector{D,Float64}` unit velocity vector
- `speed::Float64` magnitude of the velocity vector
- `motility::M` motile pattern of the microbe

and the default parameters
- `turn_rate::Float64 = 1.0` frequency of reorientations
- `rotational_diffusivity::Float64 = 0.0` coefficient of brownian rotational diffusion
- `radius::Float64 = 0.0` equivalent spherical radius of the microbe
- `state::Float64 = 0.0` generic variable for a scalar internal state
"""
@agent struct Microbe{D,M}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,M}
    speed::Float64
    motility::M
    turn_rate::Float64 = 1.0
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    state::Float64 = 0.0
end

r2dig(x) = round(x, digits=2)
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D,M}) where {D,M}
    println(io, "$(typeof(m)) with $(M)")
    println(io, "position (μm): $(r2dig.(position(m))); velocity (μm/s): $(r2dig.(velocity(m)))")
    println(io, "average unbiased turn rate (Hz): $(r2dig(turn_rate(m)))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate])
    print(io, "other properties: " * join(s, ", "))
end
