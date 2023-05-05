export neighborlist!, update_neighborlist!

"""
    neighborlist!(model::ABM, x::AbstractVector, cutoff::Real, key::Symbol)

Initialise a neighbor list between objects `x` with neighbor radius `cutoff` and
add it to the `model` properties with name `key`.
"""
function neighborlist!(model::ABM, x::AbstractVector, cutoff::Real, key::Symbol)
    model.properties[key] = neighborlist(x, model.space, cutoff)
end
"""
    neighborlist!(model::ABM, x::AbstractVector, y::AbstractVector, cutoff::Real, key::Symbol)

Initialise a neighbor list between objects `x` and `y` with neighbor radius `cutoff` and
add it to the `model` properties with name `key`.
"""
function neighborlist!(model::ABM, x::AbstractVector, y::AbstractVector, cutoff::Real, key::Symbol)
    model.properties[key] = neighborlist(x, y, model.space, cutoff)
end

neighborlist(x::AbstractVector, model::ABM, cutoff::Real) =
    neighborlist(x, model.space, cutoff)
neighborlist(x::AbstractVector, y::AbstractVector, model::ABM, cutoff::Real) =
    neighborlist(x, y, model.space, cutoff)
"""
    function neighborlist(x::AbstractVector, space::ContinuousSpace, cutoff::Real)

Initialise a neighbor list between objects in the vector `x` in `space`
using a neighbor radius `cutoff`.
"""
function neighborlist(
    x::AbstractVector, space::ContinuousSpace{D,P}, cutoff::Real
) where {D,P}
    pos_x = [SVector{length(xᵢ)}(xᵢ) for xᵢ in _pos.(x)]
    PeriodicSystem(
        xpositions = pos_x,
        # if space is not periodic the unitcell size must be extended by cutoff
        unitcell = SVector(spacesize(space) .+ (P ? 0.0 : cutoff)),
        cutoff = cutoff,
        output = 0.0
    )
end
"""
    function neighborlist(x::AbstractVector, y::AbstractVector, space::ContinuousSpace, cutoff::Real)

Initialise a neighbor list between objects in the vectors `x` and `y` in `space`
using a neighbor radius `cutoff`.
"""
function neighborlist(
    x::AbstractVector, y::AbstractVector,
    space::ContinuousSpace{D,P}, cutoff::Real
) where {D,P}
    pos_x = [SVector{length(xᵢ)}(xᵢ) for xᵢ in _pos.(x)]
    pos_y = [SVector{length(yᵢ)}(yᵢ) for yᵢ in _pos.(y)]
    PeriodicSystem(
        xpositions = pos_x,
        ypositions = pos_y,
        # if space is not periodic the unitcell size must be extended by cutoff
        unitcell = SVector(spacesize(space) .+ (P ? 0.0 : cutoff)),
        cutoff = cutoff,
        output = 0.0
    )
end


"""
    update_neighborlist!(model, listkey)

Updates the position of all microbes in the neighbor list stored in `model.properties[listkey]`.
This function assumes that the positions of microbes have been passed
as the `x` argument to `neighborlist!`.
"""
function update_neighborlist!(model, listkey)
    for microbe in allagents(model)
        update_neighborlist!(microbe, model, listkey)
    end
end
function update_neighborlist!(microbe::AbstractMicrobe, model, listkey)
    neighbor_list = model.properties[listkey]
    neighbor_list.xpositions[microbe.id] = SVector(microbe.pos)
end
