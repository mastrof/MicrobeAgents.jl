export Obstacles, Walls, obstacles, nearest_surface, random_free_position

"""
    Obstacles(space, primitives; contact = Halt(), margin = nothing, cellsize = nothing)
Static set of solid `primitives` (a vector of one concrete `AbstractObstacle` type)
in `space`, whose surfaces interact with microbes through the `contact` rule
(`Halt()`, `Slide()`, `Reflect()`, `Randomize()` or a function, see `contact!`).
Pass a `Tuple` of sets as the `:obstacles` model property.
A set is indexed and iterated as its primitives: `obstacles(model)[c.set][c.index]`
is the primitive involved in a contact `c` (or returned by `nearest_surface`).
Obstacles are ignored by `pathfinder_step!`: pathfinding and obstacles
are not meant to be used together.

Bounded primitives are indexed in a static cell list. `cellsize` (default: twice
the largest bounding radius) and `margin` (default: half the cell size) only
affect performance; queries are fastest when `margin` exceeds the microbe radius
plus the distance travelled in one step.
"""
struct Obstacles{D,P<:AbstractObstacle{D},C,I}
    primitives::Vector{P}
    contact::C
    index::I
    extent::SVector{D,Float64}
    periodic::NTuple{D,Bool}
end

_periodicity(::ContinuousSpace{D,P}) where {D,P} = P isa Bool ? ntuple(_ -> P, Val(D)) : P

function Obstacles(space::ContinuousSpace{D}, primitives::AbstractVector{<:AbstractObstacle{D}};
    contact = Halt(), margin = nothing, cellsize = nothing
) where {D}
    P = eltype(primitives)
    isempty(primitives) && throw(ArgumentError("`primitives` must not be empty"))
    isconcretetype(P) || throw(ArgumentError(
        "`primitives` must have a concrete element type, got $P; " *
        "use a separate `Obstacles` set for each primitive type"))
    prims = Vector{P}(primitives)
    extent = SVector{D,Float64}(spacesize(space))
    periodic = _periodicity(space)
    index = all(p -> isfinite(boundingradius(p)), prims) ?
        CellIndex(prims, extent, periodic, margin, cellsize) : NoIndex()
    Obstacles{D,P,typeof(contact),typeof(index)}(prims, contact, index, extent, periodic)
end

Base.getindex(set::Obstacles, i) = getindex(set.primitives, i)
Base.length(set::Obstacles) = length(set.primitives)
Base.eachindex(set::Obstacles) = eachindex(set.primitives)
Base.iterate(set::Obstacles, state...) = iterate(set.primitives, state...)
Base.eltype(::Type{<:Obstacles{D,P}}) where {D,P} = P

"""
    Walls(space; contact = Halt())
`Obstacles` made of one `ObstacleHalfSpace` per side of each non-periodic
dimension of `space`: the domain boundaries.
"""
function Walls(space::ContinuousSpace{D}; contact = Halt()) where {D}
    extent = spacesize(space)
    periodic = _periodicity(space)
    all(periodic) && throw(ArgumentError("`Walls` requires a non-periodic dimension"))
    walls = ObstacleHalfSpace{D}[]
    for d in 1:D
        periodic[d] && continue
        e = _unit(Val(D), d)
        push!(walls, ObstacleHalfSpace(zero(e), e))
        push!(walls, ObstacleHalfSpace(extent[d] * e, -e))
    end
    Obstacles(space, walls; contact)
end

# Reduce `op(acc, i, xe)` over candidate primitives `i` that may have their
# surface within `cutoff` of `x`; `xe` is the point at which to evaluate `sdf`
# (accounts for periodic images).
function _reduce_candidates(op::F, acc, set::Obstacles{D,P,C,NoIndex}, x, cutoff) where {F,D,P,C}
    # all primitives, at the periodic images of `x` (nearest images only)
    images = CartesianIndices(ntuple(d -> set.periodic[d] ? (-1:1) : (0:0), Val(D)))
    for w in images
        xw = x - SVector{D,Float64}(Tuple(w)) .* set.extent
        for i in eachindex(set.primitives)
            acc = op(acc, i, xw)
        end
    end
    return acc
end
# number of cells to visit beyond a cell boundary, for a query reaching `r` past it
_reach(r, h, n) = r ≤ 0 ? 0 : (r ≥ n * h ? n : ceil(Int, r / h))
function _reduce_candidates(op::F, acc, set::Obstacles{D,P,C,<:CellIndex}, x, cutoff) where {F,D,P,C}
    idx = set.index
    c = ntuple(d -> clamp(floor(Int, x[d] / idx.h[d]) + 1, 1, idx.dims[d]), Val(D))
    extra = cutoff - idx.margin
    # neighbor cells are visited only on the sides reached by `extra` from `x`
    klo = ntuple(d -> _reach(extra - (x[d] - (c[d] - 1) * idx.h[d]), idx.h[d], idx.dims[d]), Val(D))
    khi = ntuple(d -> _reach(extra - (c[d] * idx.h[d] - x[d]), idx.h[d], idx.dims[d]), Val(D))
    lin = LinearIndices(idx.dims)
    for o in CartesianIndices(ntuple(d -> -klo[d]:khi[d], Val(D)))
        j = c .+ Tuple(o)
        inside = true
        for d in 1:D
            set.periodic[d] || (inside &= 1 ≤ j[d] ≤ idx.dims[d])
        end
        inside || continue
        w = ntuple(d -> fld(j[d] - 1, idx.dims[d]), Val(D)) # periodic wraps
        jw = ntuple(d -> j[d] - w[d] * idx.dims[d], Val(D))
        xw = x - SVector{D,Float64}(w) .* set.extent
        l = lin[jw...]
        for e in idx.start[l]:idx.start[l+1]-1
            i, s = idx.entries[e]
            acc = op(acc, i, xw - s)
        end
    end
    return acc
end

"""
    nearest_in_set(set::Obstacles, x, cutoff) -> (d, n, i)
Signed distance `d`, outward normal `n` and index `i` of the primitive of `set`
with the closest surface within `cutoff` of `x`; `(Inf, 0, 0)` if none.
"""
function nearest_in_set(set::Obstacles{D}, x::SVector{D}, cutoff::Real) where {D}
    prims = set.primitives
    # closure locals must not shadow outer locals (would box them and allocate)
    op = (acc, j, xe) -> (dj = sdf(prims[j], xe); dj < acc[1] ? (dj, j, xe) : acc)
    dmin, i, xmin = _reduce_candidates(op, (Inf, 0, x), set, x, cutoff)
    (i == 0 || dmin > cutoff) && return (Inf, zero(SVector{D,Float64}), 0)
    return (dmin, surface_normal(prims[i], xmin), i)
end

# nearest over a tuple of sets: (d, n, set index, primitive index)
@inline _nearest(::Tuple{}, k, x, cutoff, best) = best
@inline function _nearest(sets::Tuple, k, x, cutoff, best)
    d, n, i = nearest_in_set(first(sets), x, cutoff)
    best = d < best[1] ? (d, n, k, i) : best
    _nearest(Base.tail(sets), k + 1, x, cutoff, best)
end
_nearest(sets::Tuple, x::SVector{D}, cutoff) where {D} =
    _nearest(sets, 1, x, cutoff, (Inf, zero(SVector{D,Float64}), 0, 0))

"""
    obstacles(model)
Returns the `Tuple` of `Obstacles` sets of `model`.
"""
obstacles(model::ABM) = abmproperties(model).obstacles

"""
    nearest_surface(microbe, model; cutoff = Inf) -> (; distance, normal, set, index)
    nearest_surface(x, model; cutoff = Inf)
Closest obstacle surface to `microbe` (or to point `x`) within `cutoff`: `distance`
from the microbe surface (center distance minus `radius(microbe)`), outward `normal`,
index of the obstacle `set` in `obstacles(model)` and of the primitive (`index`) in the set.
`distance = Inf` (and indices `0`) if there is none.
In per-step behaviors pass a finite `cutoff`: with `cutoff = Inf` the query may
visit every cell of the obstacle index.
"""
function nearest_surface(microbe::AbstractMicrobe, model::ABM; cutoff::Real = Inf)
    R = radius(microbe)
    q = nearest_surface(position(microbe), model; cutoff = R + cutoff)
    return (; q..., distance = q.distance - R)
end
function nearest_surface(x::Union{SVector{D},NTuple{D,Real}}, model::ABM{<:ContinuousSpace{D}};
    cutoff::Real = Inf
) where {D}
    d, n, k, i = _nearest(obstacles(model), SVector{D,Float64}(Tuple(x)), cutoff)
    return (distance = d, normal = n, set = k, index = i)
end

"""
    random_free_position(model; radius = 0.0, maxtries = 10_000)
Random position in the domain at distance `≥ radius` from all obstacle surfaces
(outside all obstacles). Throws an error after `maxtries` rejected samples.
"""
function random_free_position(model::ABM{<:ContinuousSpace{D}};
    radius::Real = 0.0, maxtries::Int = 10_000
) where {D}
    for _ in 1:maxtries
        x = SVector{D,Float64}(random_position(model))
        d, = _nearest(obstacles(model), x, radius)
        d ≥ radius && return x
    end
    error("no free position found after $maxtries attempts")
end

# validation of the `:obstacles` model property
_check_obstacles(o::Obstacles, space) = _check_obstacles((o,), space)
function _check_obstacles(obs::Tuple, space::ContinuousSpace{D}) where {D}
    for o in obs
        o isa Obstacles{D} || throw(ArgumentError(
            "`obstacles` must be a Tuple of `Obstacles` with the dimensionality of the space"))
        (o.extent == SVector{D,Float64}(spacesize(space)) && o.periodic == _periodicity(space)) ||
            throw(ArgumentError("`Obstacles` were built for a different space (extent or periodicity)"))
    end
    return obs
end
_check_obstacles(o, space) = throw(ArgumentError(
    "`obstacles` must be a Tuple of `Obstacles`, got $(typeof(o))"))
