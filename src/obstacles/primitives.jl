export AbstractObstacle, ObstacleSphere, ObstacleHalfSpace, ObstacleBox, ObstacleCapsule, ObstacleSDF
export sdf

"""
    AbstractObstacle{D}
Supertype of solid primitives in `D` dimensions, collected into `Obstacles`.

Required method: `sdf(p, x)`.
Optional, to enable spatial indexing:
`MicrobeAgents.center(p)` and `MicrobeAgents.boundingradius(p)`
(a ball enclosing the solid; default radius `Inf`, i.e. not indexed).
Optional:
`MicrobeAgents.surface_normal(p, x)` (default: finite differences of `sdf`) and
`MicrobeAgents.isconvex(::Type{<:MyObstacle}) = true` for convex solids
(exact collisions; non-convex solids are resolved approximately).
"""
abstract type AbstractObstacle{D} end

"""
    sdf(p::AbstractObstacle, x)
Signed distance from point `x` to the surface of `p` (negative inside the solid).
"""
function sdf end

"""
    surface_normal(p::AbstractObstacle, x)
Outward unit normal of the surface of `p` nearest to `x` (gradient of `sdf`).
Default: central finite differences of `sdf`.
"""
function surface_normal(p::AbstractObstacle{D}, x::SVector{D}) where {D}
    h = cbrt(eps()) * max(1.0, maximum(abs, x)) # optimal finite difference step
    g = SVector{D,Float64}(ntuple(Val(D)) do i
        e = SVector{D,Float64}(ntuple(j -> j == i ? h : 0.0, Val(D)))
        (sdf(p, x + e) - sdf(p, x - e)) / 2h
    end)
    return g / norm(g)
end

center(p::AbstractObstacle{D}) where {D} = zero(SVector{D,Float64})
boundingradius(p::AbstractObstacle) = Inf
isconvex(::Type{<:AbstractObstacle}) = false
isconvex(p::AbstractObstacle) = isconvex(typeof(p))

_svec(D, x) = SVector{D,Float64}(Tuple(x))
# unit vector along the first axis, fallback for undefined normals
_unit(::Val{D}, d = 1) where {D} = SVector{D,Float64}(ntuple(i -> i == d ? 1.0 : 0.0, Val(D)))

"""
    ObstacleSphere(center, radius)
Solid sphere (disk in 2D, segment in 1D).
"""
struct ObstacleSphere{D} <: AbstractObstacle{D}
    center::SVector{D,Float64}
    radius::Float64
end
ObstacleSphere(center, radius::Real) =
    ObstacleSphere{length(center)}(_svec(length(center), center), Float64(radius))
sdf(p::ObstacleSphere{D}, x::SVector{D}) where {D} = norm(x - p.center) - p.radius
function surface_normal(p::ObstacleSphere{D}, x::SVector{D}) where {D}
    r = x - p.center
    l = norm(r)
    iszero(l) ? _unit(Val(D)) : r / l
end
center(p::ObstacleSphere) = p.center
boundingradius(p::ObstacleSphere) = p.radius
isconvex(::Type{<:ObstacleSphere}) = true

"""
    ObstacleHalfSpace(point, normal)
Solid half-space bounded by the plane through `point`;
`normal` points away from the solid.
Unbounded, hence never spatially indexed.
"""
struct ObstacleHalfSpace{D} <: AbstractObstacle{D}
    point::SVector{D,Float64}
    normal::SVector{D,Float64}
    # inner constructor: no implicit `ObstacleHalfSpace(::SVector, ::SVector)`
    # method bypassing the normalization below
    ObstacleHalfSpace{D}(point, normal) where {D} = new{D}(point, normal)
end
function ObstacleHalfSpace(point, normal)
    D = length(point)
    n = _svec(D, normal)
    ObstacleHalfSpace{D}(_svec(D, point), n / norm(n))
end
sdf(p::ObstacleHalfSpace{D}, x::SVector{D}) where {D} = dot(x - p.point, p.normal)
surface_normal(p::ObstacleHalfSpace{D}, x::SVector{D}) where {D} = p.normal
isconvex(::Type{<:ObstacleHalfSpace}) = true

"""
    ObstacleBox(lo, hi)
Solid axis-aligned box with opposite corners `lo` and `hi`.
"""
struct ObstacleBox{D} <: AbstractObstacle{D}
    center::SVector{D,Float64}
    halfwidth::SVector{D,Float64}
    # inner constructor: no implicit `ObstacleBox(center, halfwidth)` method
    # shadowing `ObstacleBox(lo, hi)` for `SVector` arguments
    ObstacleBox{D}(center, halfwidth) where {D} = new{D}(center, halfwidth)
end
function ObstacleBox(lo, hi)
    D = length(lo)
    l, h = _svec(D, lo), _svec(D, hi)
    all(l .< h) || throw(ArgumentError("`lo` must be smaller than `hi` along all dimensions"))
    ObstacleBox{D}((l + h) / 2, (h - l) / 2)
end
function sdf(p::ObstacleBox{D}, x::SVector{D}) where {D}
    q = abs.(x - p.center) - p.halfwidth
    norm(max.(q, 0.0)) + min(maximum(q), 0.0)
end
function surface_normal(p::ObstacleBox{D}, x::SVector{D}) where {D}
    r = x - p.center
    q = abs.(r) - p.halfwidth
    σ = map(v -> v < 0 ? -1.0 : 1.0, r)
    if maximum(q) > 0 # outside: direction from the closest point of the box
        g = max.(q, 0.0) .* σ
        return g / norm(g)
    else # inside: the closest face
        k = argmax(q)
        return SVector{D,Float64}(ntuple(i -> i == k ? σ[i] : 0.0, Val(D)))
    end
end
center(p::ObstacleBox) = p.center
boundingradius(p::ObstacleBox) = norm(p.halfwidth)
isconvex(::Type{<:ObstacleBox}) = true

"""
    ObstacleCapsule(a, b, radius)
Solid capsule: all points within `radius` of the segment from `a` to `b`
(rods, fibres).
"""
struct ObstacleCapsule{D} <: AbstractObstacle{D}
    a::SVector{D,Float64}
    b::SVector{D,Float64}
    radius::Float64
end
ObstacleCapsule(a, b, radius::Real) = ObstacleCapsule{length(a)}(
    _svec(length(a), a), _svec(length(a), b), Float64(radius)
)
function _capsule_offset(p::ObstacleCapsule{D}, x::SVector{D}) where {D}
    ba = p.b - p.a
    l2 = dot(ba, ba)
    t = iszero(l2) ? 0.0 : clamp(dot(x - p.a, ba) / l2, 0.0, 1.0)
    x - (p.a + t * ba)
end
sdf(p::ObstacleCapsule{D}, x::SVector{D}) where {D} = norm(_capsule_offset(p, x)) - p.radius
function surface_normal(p::ObstacleCapsule{D}, x::SVector{D}) where {D}
    v = _capsule_offset(p, x)
    l = norm(v)
    iszero(l) ? _unit(Val(D)) : v / l
end
center(p::ObstacleCapsule) = (p.a + p.b) / 2
boundingradius(p::ObstacleCapsule) = norm(p.b - p.a) / 2 + p.radius
isconvex(::Type{<:ObstacleCapsule}) = true

"""
    ObstacleSDF(f, center, boundingradius = Inf; convex = false)
Solid defined by an arbitrary signed distance function `f(x::SVector{D}) -> Real`
(negative inside; it must not overestimate the distance to the surface).
`center` and `boundingradius` define a ball enclosing the solid, used for spatial
indexing (`Inf`: not indexed). Set `convex = true` only if the solid is convex:
collisions with convex solids are exact, while non-convex solids are found by ray
marching, which may miss features thinner than 1% of the larger of the microbe radius
and the step length, and leaves shallow overlaps (removed at the end of the step)
when sliding along concave parts (prefer unions of convex primitives).
Normals are computed by finite differences.
"""
struct ObstacleSDF{D,F,C} <: AbstractObstacle{D}
    f::F
    center::SVector{D,Float64}
    boundingradius::Float64
end
ObstacleSDF(f, center, boundingradius::Real = Inf; convex::Bool = false) =
    ObstacleSDF{length(center),typeof(f),convex}(
        f, _svec(length(center), center), Float64(boundingradius)
    )
sdf(p::ObstacleSDF{D}, x::SVector{D}) where {D} = Float64(p.f(x))
center(p::ObstacleSDF) = p.center
boundingradius(p::ObstacleSDF) = p.boundingradius
isconvex(::Type{<:ObstacleSDF{D,F,C}}) where {D,F,C} = C
