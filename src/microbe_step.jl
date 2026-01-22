export microbe_step!, microbe_pathfinder_step!
export move_step!, move_step_pathfinder!, affect_step!, reorient_step!
export switching_probability, can_turn

"""
    microbe_step!(microbe, model)
Perform an integration step for `microbe`. In order:
1. Update `microbe` position according to its current velocity
2. Reorient `microbe` through rotational diffusion
3. Update internal state of `microbe` through a customizable `affect!` function
4. Reorient `microbe` (if the motile state has a reorientation component)
5. Update motile state of microbe (transitions occur as Poisson events)

Step 1 and 2 are implemented by `move_step!`.
Step 3 is implemented by `affect_step!`.
Step 4 and 5 are implemented by `reorient_step!`.
"""
function microbe_step!(microbe::AbstractMicrobe, model::ABM)
    move_step!(microbe, model)
    affect_step!(microbe, model)
    reorient_step!(microbe, model)
end

"""
    microbe_pathfinder_step!(microbe, model)
Identical to `microbe_step!` but motion is constrained by the pathfinder
(`model.pathfinder`) which defines inaccessible regions of space.
"""
function microbe_pathfinder_step!(microbe::AbstractMicrobe, model::ABM)
    move_step_pathfinder!(microbe, model)
    affect_step!(microbe, model)
    reorient_step!(microbe, model)
end

"""
    move_agent!(microbe::AbstractMicrobe, model::ABM{<:ContinuousSpace}, dt::Real)
Propagate the microbe forward one step according to its velocity vector
(after updating its velocity with `update_vel!` if configured).

This function respects the space size.
"""
function Agents.move_agent!(microbe::AbstractMicrobe, model::ABM{<:ContinuousSpace}, dt::Real)
    abmspace(model).update_vel!(microbe, model)
    walk!(microbe, velocity(microbe).*dt, model)
end

"""
    move_step!(microbe::AbstractMicrobe, model::ABM)
Subroutine to perform one step of translation followed by
reorientation due to rotational diffusion.
"""
function move_step!(microbe::AbstractMicrobe, model::ABM)
    dt = abmtimestep(model)
    move_agent!(microbe, model, dt) # translation
    rotational_diffusion!(microbe, model) # angular noise
end

"""
    move_step_pathfinder!(microbe::AbstractMicrobe, model::ABM)
Subroutine to perform one step of translation constrained
by `model.pathfinder`, followed by reorientation due to
rotational diffusion.
"""
function move_step_pathfinder!(microbe::AbstractMicrobe, model::ABM)
    dt = abmtimestep(model)
    pathfinder_step!(microbe, model, dt) # translation
    rotational_diffusion!(microbe, model) # angular noise
end

"""
    affect_step!(microbe::AbstractMicrobe, model::ABM)
Subroutine to update the internal state of the microbe.
"""
function affect_step!(microbe::AbstractMicrobe, model::ABM)
    model.affect!(microbe, model)
end

"""
    reorient_step!(microbe::AbstractMicrobe, model::ABM)
Subroutine to perform active reorientation and/or
switch the microbe to a new motile state.
"""
function reorient_step!(microbe::AbstractMicrobe, model::ABM)
    if can_turn(microbe)
        turn!(microbe, model)
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
    switching_probability(microbe::AbstractMicrobe, model::ABM)
Evaluate the probability of switching to another motile state.

The probability is defined as β*dt/τ, where:
- τ is the average duration of current motile state
- dt is the integration timestep (`abmtimestep(model)`)
- β is a physiological switching "bias" (`bias(microbe)`).
  This is defined differently for each microbe type.
  An unbiased state corresponds to β=1.
  β>1 implies higher probability to change motile state.
"""
function switching_probability(microbe::AbstractMicrobe, model::ABM)
    dt = abmtimestep(model)
    M = motilestate(microbe)
    τ = duration(M)
    β = biased(M) ? bias(microbe) : 1.0
    return β * dt / τ
end

"""
    can_turn(microbe::AbstractMicrobe)
Check whether the microbe is in a `MotileState`
associated with a non-null reorientation angle.
"""
function can_turn(microbe::AbstractMicrobe)
    m = motilestate(microbe)
    !isnothing(angle(m))
end
