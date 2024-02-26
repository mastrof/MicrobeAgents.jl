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
function microbe_step!(microbe::AbstractMicrobe, model)
    dt = model.timestep # integration timestep
    # update microbe position
    move_agent!(microbe, model, speed(microbe)*dt)
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
