export Arccos

struct Arccos{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    Arccos{T}(a::Real, b::Real) where {T} = new{T}(T(a), T(b))
end
function Arccos(a::Real, b::Real; check_args::Bool=true)
    Distributions.@check_args Arccos a<b -1≤a≤1 -1≤b≤1
    return Arccos{Float64}(Float64(a), Float64(b))
end # function
Arccos() = Arccos(-1,1)
Base.rand(rng::AbstractRNG, d::Arccos) = acos(rand(rng, Uniform(d.a, d.b)))