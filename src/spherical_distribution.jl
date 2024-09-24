export Spherical

struct Spherical{T<:Real} <: ContinuousUnivariateDistribution
    θmin::T
    θmax::T
end
"""
    Spherical(θmin, θmax)
Uniform distribution of angles on the spherical sector
with aperture defined by the angles `θmin` and `θmax` (in radians).

The default values, `θmin=0` and `θmax=π`, give the uniform
distribution over the entire sphere.
"""
function Spherical(θmin::Real=0, θmax::Real=π)
    Spherical{Float64}(Float64(θmin), Float64(θmax))
end
function Base.rand(rng::AbstractRNG, d::Spherical)
    u = rand(rng, Uniform(-1, +1))
    d.θmin + ((d.θmax-d.θmin) / π) * acos(u)
end

