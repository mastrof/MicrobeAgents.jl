export Microbe, chemotaxis!, bias

"""
    Microbe{D,N} <: AbstractMicrobe{D,N}
Base microbe type for simple simulations.
`D` is the space dimensionality and `N` is the number of
states of the microbe's motility pattern.

`Microbe` has the required fields
- `id::Int` an identifier used internally
- `pos::SVector{D,Float64}` spatial position
- `vel::SVector{D,Float64}` unit velocity vector
- `speed::Float64` magnitude of the velocity vector
- `motility::Motility{N}` motility pattern

and the default parameters
- `rotational_diffusivity::Float64 = 0.0` coefficient of brownian rotational diffusion
- `radius::Float64 = 0.0` equivalent spherical radius of the microbe
- `state::Float64 = 0.0` generic variable for a scalar internal state
"""
@agent struct Microbe{D,N}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N}
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    state::Float64 = 0.0
end

# fallback functions for default random behavior
"""
    chemotaxis!(microbe::AbstractMicrobe)
Modify the internal state of the microbe using chemotaxis.
Each microbe type with chemotactic capabilities has a distinct
implementation of this function, corresponding to different
models of chemotaxis.
For the base `Microbe` type, this does nothing.
"""
chemotaxis!(microbe::AbstractMicrobe, model) = nothing
"""
    bias(microbe::AbstractMicrobe)
Return the current bias in the switching probability of the microbe
(see `switching_probability`).
For the base, non-chemotactic `Microbe` type, this is just 1.
For other microbe types, the bias is evaluated based on the
specifics of the chemotactic model.
A bias>1 indicates higher likelihood to switch motile state.
"""
bias(microbe::AbstractMicrobe) = 1.0


r2dig(x) = round(x, digits=2)
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D,N}) where {D,N}
    println(io, "$(typeof(m)) with $(N)-state motility pattern")
    println(io, "position (μm): $(r2dig.(position(m))); velocity (μm/s): $(r2dig.(velocity(m)))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate])
    print(io, "other properties: " * join(s, ", "))
end
