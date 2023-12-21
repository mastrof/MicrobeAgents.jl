export random_velocity, random_speed, ChainedFunction, →

"""
    random_velocity(model)
Generate a random velocity vector with unit norm respecting the
dimensionality of `model`.
"""
function random_velocity(model::AgentBasedModel{<:ContinuousSpace{D}}) where D
    random_velocity(abmrng(model), D)
end
function random_velocity(rng::AbstractRNG, D::Int)
    v = rand(rng, SVector{D}) .* 2 .- 1
    mag = sqrt(dot(v,v))
    return v ./ mag
end

"""
    random_speed(microbe, model)
Generate a random speed from the motile pattern of `microbe`.
"""
function random_speed(microbe::AbstractMicrobe, model::AgentBasedModel)
    rand(abmrng(model), speed(motilepattern(microbe)))
end

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
