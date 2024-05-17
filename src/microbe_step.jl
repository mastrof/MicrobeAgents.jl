export microbe_step!, microbe_pathfinder_step!

"""
    microbe_step!(microbe, model)
Perform an integration step for `microbe`. In order:
1. Update `microbe` position according to its current velocity
2. Reorient `microbe` through rotational diffusion
3. Update the state of `microbe` through a dedicated `affect!` function
4. Evaluate instantaneous turn rate through dedicated `turnrate` function
5. Reorient `microbe` if necessary
"""
function microbe_step!(microbe::AbstractMicrobe, model::AgentBasedModel)
    dt = abmtimestep(model)
    move_agent!(microbe, model, speed(microbe)*dt) # translation
    rotational_diffusion!(microbe, model) # rotational noise
    turn!(microbe, model) # reorientation
    model.affect!(microbe, model) # update microbe's internal state
    # probability to switch motile state
    p = dt * bias(microbe) / duration(motilepattern(microbe))
    if rand(abmrng(model)) < p
        update_motilestate!(microbe, model)
    end
end

"""
    microbe_pathfinder_step!(microbe, model)
Identical to `microbe_step!` except that motion is constrained by the pathfinder
(`model.pathfinder`) which defines inaccessible regions of space.
"""
function microbe_pathfinder_step!(microbe::AbstractMicrobe, model)
    dt = model.timestep # integration timestep
    # update microbe position
    pathfinder_step!(microbe, model, dt)
    # reorient through rotational diffusion
    rotational_diffusion!(microbe, model)
    # update microbe state
    model.affect!(microbe, model)
    # evaluate instantaneous turn rate
    ω = turnrate(microbe) * tumblebias(microbe)
    if rand(abmrng(model)) < ω * dt # if true reorient microbe
        turn!(microbe, model)
    end
    nothing
end
