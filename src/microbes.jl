export Microbe

"""
    Microbe{D,N,B} <: AbstractMicrobe{D,N}
Microbe type of MicrobeAgents.jl.
`D` is the space dimensionality, `N` the number of motile states,
`B` the type of the `behaviors` container (inferred, never written by hand).

Fields:
- `id::Int` identifier used internally
- `pos::SVector{D,Float64}` position
- `vel::SVector{D,Float64}` unit direction of motion
- `speed::Float64` sampled speed of the current motile state (see `speed`)
- `motility::Motility{N}` motility pattern
- `rotational_diffusivity::Float64 = 0.0` rotational diffusion coefficient
- `radius::Float64 = 0.0` equivalent spherical radius
- `behaviors::B = ()` `Tuple` or `NamedTuple` of behaviors (see `Behavior`)
"""
@agent struct Microbe{D,N,B}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N}
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    behaviors::B = ()
end
Microbe{D,N}(; behaviors = (), kwargs...) where {D,N} =
    Microbe{D,N,typeof(behaviors)}(; behaviors, kwargs...)

behaviors(m::Microbe) = m.behaviors

r2dig(x) = round(x, digits=2)
_label(b) = string(nameof(typeof(b)))
_label(f::Function) = string(nameof(f))
_labels(bs::Tuple) = map(_label, bs)
_labels(bs::NamedTuple) = map((k, b) -> "$k = $(_label(b))", keys(bs), values(bs))
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D,N}) where {D,N}
    println(io, "$(nameof(typeof(m))){$D} with $(N)-state motility pattern")
    println(io, "position (μm): $(r2dig.(position(m))); velocity (μm/s): $(r2dig.(velocity(m)))")
    bs = behaviors(m)
    print(io, "behaviors: ", isempty(bs) ? "none" : join(_labels(bs), ", "))
end
