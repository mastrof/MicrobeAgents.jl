export pathfinder!, initialise_pathfinder, pathfinder_step!

function pathfinder!(model::AgentBasedModel, walkmap::BitArray)
    pathfinder = initialise_pathfinder(abmspace(model), walkmap)
    abmproperties(model)[:pathfinder] = pathfinder
end

function initialise_pathfinder(
    extent::Real, periodic::Bool,
    walkmap::BitArray{D}
) where D
    initialise_pathfinder(ntuple(_->extent,D), periodic, walkmap)
end
function initialise_pathfinder(
    extent::NTuple{D,<:Real}, periodic::Bool,
    walkmap::BitArray{D}
) where D
    space = ContinuousSpace(extent; periodic)
    AStar(space; walkmap)
end
function initialise_pathfinder(space::ContinuousSpace{D}, walkmap::BitArray{D}) where D
    AStar(space; walkmap)
end

"""
    pathfinder_step!(microbe::AbstractMicrobe, model::AgentBasedModel, dt::Real)
Perform an integration step for `microbe` motion with pathfinding
(`model.pathfinder`).
"""
function pathfinder_step!(microbe::AbstractMicrobe, model::AgentBasedModel, dt::Real)
    U = microbe.speed
    target_position = @. microbe.pos + U*microbe.vel*dt
    plan_route!(microbe, target_position, model.pathfinder)
    move_along_route!(microbe, model, model.pathfinder, U, dt)
end
