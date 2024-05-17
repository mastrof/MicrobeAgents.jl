export microbe_step!, microbe_pathfinder_step!

"""
    microbe_step!(microbe, model)
Perform an integration step for `microbe`. In order:
1. Update `microbe` position according to its current velocity
2. Reorient `microbe` through rotational diffusion
3. Update internal state of `microbe` through a customizable `affect!` function
4. Reorient `microbe` (if the motile state has a reorientation component)
5. Update motile state of microbe (transitions occur as Poisson events)
"""
function microbe_step!(microbe::AbstractMicrobe, model::AgentBasedModel)
    dt = abmtimestep(model)
    move_agent!(microbe, model, speed(microbe)*dt) # translation
    rotational_diffusion!(microbe, model) # rotational noise
    model.affect!(microbe, model) # update microbe's internal state
    turn!(microbe, model) # reorientation
    p = switching_probability(microbe, model)
    if rand(abmrng(model)) < p
        # transition to new motile state
        update_motilestate!(microbe, model)
    end
end

function switching_probability(microbe, model)
    dt = abmtimestep(model)
    M = motilestate(microbe)
    τ = duration(M)
    β = biased(M) ? bias(microbe) : 1.0
    return β * dt / τ
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
