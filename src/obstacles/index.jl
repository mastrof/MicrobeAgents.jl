# Static spatial index of obstacle primitives (obstacle-side cell list).

struct NoIndex end

"""
    CellIndex{D}
Uniform grid over the domain (CSR storage). Cell `j` lists `(i, s)` entries:
primitive `i`, translated by the periodic image shift `s`, may have its surface
within `margin` of the cell.
"""
struct CellIndex{D}
    h::SVector{D,Float64}
    dims::NTuple{D,Int}
    margin::Float64
    start::Vector{Int}
    entries::Vector{Tuple{Int,SVector{D,Float64}}}
end

const MAXCELLS = 1_000_000

function CellIndex(prims::AbstractVector{<:AbstractObstacle{D}},
    extent::SVector{D,Float64}, periodic::NTuple{D,Bool}, margin, cellsize
) where {D}
    rmax = maximum(boundingradius, prims)
    h0 = max(isnothing(cellsize) ? 2rmax : Float64(cellsize), (prod(extent) / MAXCELLS)^(1 / D))
    # integer number of cells per dimension, so that periodic wraps align with cells
    dims = ntuple(d -> max(1, floor(Int, extent[d] / h0)), Val(D))
    h = SVector{D,Float64}(ntuple(d -> extent[d] / dims[d], Val(D)))
    m = isnothing(margin) ? minimum(h) / 2 : Float64(margin)
    m ≥ 0 || throw(ArgumentError("`margin` must be non-negative"))
    lists = [Tuple{Int,SVector{D,Float64}}[] for _ in CartesianIndices(dims)]
    shifts = Iterators.product(ntuple(d -> periodic[d] ? (-extent[d], 0.0, extent[d]) : (0.0,), Val(D))...)
    for (i, p) in pairs(prims), s in shifts
        sv = SVector{D,Float64}(s)
        c = center(p) + sv
        rb = boundingradius(p)
        r = rb + m
        lo = ntuple(d -> clamp(floor(Int, (c[d] - r) / h[d]) + 1, 1, dims[d]), Val(D))
        hi = ntuple(d -> clamp(floor(Int, (c[d] + r) / h[d]) + 1, 1, dims[d]), Val(D))
        for J in CartesianIndices(ntuple(d -> lo[d]:hi[d], Val(D)))
            blo = (SVector{D,Float64}(Tuple(J)) .- 1) .* h
            q = clamp.(c, blo, blo + h) # closest point of the cell to the center
            # the bounding ball contains the solid: conservative for any shape,
            # and keeps interior cells (queries from inside solids stay correct)
            norm(q - c) - rb ≤ m && push!(lists[J], (i, sv))
        end
    end
    start = Vector{Int}(undef, length(lists) + 1)
    start[1] = 1
    for l in eachindex(lists)
        start[l+1] = start[l] + length(lists[l])
    end
    CellIndex{D}(h, dims, m, start, reduce(vcat, vec(lists)))
end
