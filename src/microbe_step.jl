export microbe_step!, microbe_pathfinder_step!, switching_probability

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
    if can_turn(microbe)
        turn!(microbe, model) # reorientation
    end
    p = switching_probability(microbe, model)
    if rand(abmrng(model)) < p
        # transition to new motile state
        update_motilestate!(microbe, model)
        # if the new motile state is a TurnState with duration 0
        # immediately turn and transition to another state
        new_motilestate = motilestate(microbe)
        if variantof(new_motilestate) === TurnState && iszero(duration(new_motilestate))
            turn!(microbe, model)
            update_motilestate!(microbe, model)
        end
        # sample new speed
        update_speed!(microbe, model)
    end
end

"""
    microbe_pathfinder_step!(microbe, model)
Identical to `microbe_step!` except that motion is constrained by the pathfinder
(`model.pathfinder`) which defines inaccessible regions of space.
"""
function microbe_pathfinder_step!(microbe::AbstractMicrobe, model::AgentBasedModel)
    dt = abmtimestep(model)
    move_agent!(microbe, model, dt) # translation
    rotational_diffusion!(microbe, model) # rotational noise
    model.affect!(microbe, model) # update microbe's internal state
    if can_turn(microbe)
        turn!(microbe, model) # reorientation
    end
    p = switching_probability(microbe, model)
    if rand(abmrng(model)) < p
        # transition to new motile state
        update_motilestate!(microbe, model)
        # if the new motile state is a TurnState with duration 0
        # immediately turn and transition to another state
        new_motilestate = motilestate(microbe)
        if variantof(new_motilestate) === TurnState && iszero(duration(new_motilestate))
            turn!(microbe, model)
            update_motilestate!(microbe, model)
        end
        # sample new speed
        update_speed!(microbe, model)
    end
end

function switching_probability(microbe, model)
    dt = abmtimestep(model)
    M = motilestate(microbe)
    τ = duration(M)
    β = biased(M) ? bias(microbe) : 1.0
    return β * dt / τ
end

function can_turn(microbe)
    m = motilestate(microbe)
    !isnothing(polar(m))
end
