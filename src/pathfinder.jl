export AStar # re-export from Agents.Pathfinding
export pathfinder_step!

"""
    pathfinder_step!(microbe::AbstractMicrobe, model::ABM, dt::Real)
Perform an integration step for `microbe` motion with pathfinding
(`model.pathfinder`).
"""
function pathfinder_step!(microbe::AbstractMicrobe, model::ABM, dt::Real)
    U = speed(microbe)
    target_position = normalize_position(
        position(microbe) .+ velocity(microbe).*dt,
        model
    )
    plan_route!(microbe, target_position, model.pathfinder)
    move_along_route!(microbe, model, model.pathfinder, U, dt)
end
