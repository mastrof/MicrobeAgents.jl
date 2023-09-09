export rand_vel, ChainedFunction, →, distance, distancevector

"""
    random_velocity(model)
Generate a random velocity vector with unit norm respecting the
dimensionality of `model`.
"""
function random_velocity(model::AgentBasedModel{S,A}) where {S,D,A<:AbstractMicrobe{D}}
    v = rand(abmrng(model), SVector{D}) .* 2 .- 1
    mag = sqrt(dot(v,v))
    return v ./ mag
end

"""
    random_speed(microbe, model)
Generate a random speed from the motile pattern of `microbe`.
"""
function random_speed(microbe::AbstractMicrobe, model::AgentBasedModel)
    random_speed(abmrng(model), motilepattern(microbe))
end

motilepattern(microbe::AbstractMicrobe) = microbe.motility


struct ChainedFunction{H,T} <: Function
    head::H
    tail::T
    ChainedFunction{H,T}(head, tail) where {H,T} = new{H,T}(head, tail)
    ChainedFunction(head, tail) = new{typeof(head),typeof(tail)}(head, tail)
end
(c::ChainedFunction)(x...; kw...) = (c.head(x...; kw...); c.tail(x...; kw...))
→(f) = f
→(f,g) = ChainedFunction(f,g)
→(f, g, h...) = →(f→g, h...)

function Base.show(io::IO, ::MIME"text/plain", f::ChainedFunction)
    s = join(funlist(f), " → ")
    print(io, "ChainedFunction: $s")
end
funlist(f::ChainedFunction{<:ChainedFunction,<:ChainedFunction}) = (funlist(f.head)..., funlist(f.tail)...)
funlist(f::ChainedFunction{<:Function,<:ChainedFunction}) = (f.head, funlist(f.tail)...)
funlist(f::ChainedFunction{<:Function,<:Function}) = (f.head, f.tail)

@inline _pos(a::AbstractMicrobe) = a.pos
@inline _pos(a::NTuple{D}) where D = a
"""
    distance(a, b, model)
Evaluate the euclidean distance between `a` and `b` respecting the spatial
properties of `model`.
"""
@inline distance(a, b, model) = euclidean_distance(_pos(a), _pos(b), model)
"""
    distancevector(a, b, model)
Evaluate the distance vector from `a` to `b` respecting the spatial
properties of `model`.
"""
@inline distancevector(a, b, model) = distancevector(_pos(a), _pos(b), model)
@inline function distancevector(a::NTuple{D}, b::NTuple{D}, model) where D
    extent = spacesize(model)
    ntuple(i -> wrapcoord(a[i], b[i], extent[i]), D)
end
function wrapcoord(x₁, x₂, d)
    α = (x₂-x₁)/d
    (α-round(α))*d
end
