export neighborlist!, update_neighborlist!

"""
    neighborlist!(model::ABM, cutoff::Real, key::Symbol)

Initialise a neighbor list of the microbes in `model` with neighbor radius `cutoff` and
add it to the `model` properties with name `key`.

*Neighbor lists are not supported in 1-dimensional models.*
"""
function neighborlist!(model::ABM, cutoff::Real, key::Symbol)
    D = length(spacesize(model))
    if D == 1
        throw(ArgumentError("Neighbor lists are not supported in 1-dimensional models."))
    end
    abmproperties(model)[key] = neighborlist(model, cutoff)
end

"""
    neighborlist!(model::ABM, y, cutoff::Real, key::Symbol)

Initialise a neighbor list between the microbes in `model` and the objects `y`
with neighbor radius `cutoff` and add it to the `model` properties with name `key`.

*Neighbor lists are not supported in 1-dimensional models.*
"""
function neighborlist!(model::ABM, y, cutoff, key)
    D = length(spacesize(model))
    if D == 1
        throw(ArgumentError("Neighbor lists are not supported in 1-dimensional models."))
    end
    abmproperties(model)[key] = neighborlist(model, y, cutoff)
end


function neighborlist(model::ABM, cutoff)
    microbes = make_position_vector(model)
    neighborlist(microbes, model.space, cutoff)
end
function neighborlist(model::ABM, y, cutoff)
    microbes = make_position_vector(model)
    neighborlist(microbes, make_position_vector(y), cutoff)
end


"""
    function neighborlist(x::AbstractVector, space::ContinuousSpace, cutoff::Real)

Initialise a neighbor list between objects in the vector `x` in `space`
using a neighbor radius `cutoff`.
"""
function neighborlist(
    x::AbstractVector, space::ContinuousSpace{D,P}, cutoff::Real
) where {D,P}
    PeriodicSystem(
        xpositions = x,
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
    PeriodicSystem(
        xpositions = x,
        ypositions = y,
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
    neighbor_list = abmproperties(model)[listkey]
    neighbor_list.xpositions[microbe.id] = SVector(microbe.pos)
end


@inline function make_position_vector(model::UnremovableABM)
    [SVector(_pos(a)) for a in allagents(model)]
end
@inline function make_position_vector(model::StandardABM)
    ids = sort(collect(allids(model)))
    [SVector(_pos(model[i])) for i in ids]
end
@inline function make_position_vector(x::AbstractVector)
    [SVector{length(xᵢ)}(_pos(xᵢ)) for xᵢ in x]
end
@inline function make_position_vector(x::Dict)
    indices = sort(collect(keys(x)))
    [SVector{length(x[i])}(_pos(x[i])) for i in indices]
end
