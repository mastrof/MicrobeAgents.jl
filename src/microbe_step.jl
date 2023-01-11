export microbe_step!, turnrate, affect!

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
    move_agent!(microbe, model, dt)
    # reorient through rotational diffusion
    rotational_diffusion!(microbe, model)
    # update microbe state
    affect!(microbe, model)
    # evaluate instantaneous turn rate
    ω = turnrate(microbe, model)
    if rand(model.rng) < ω*dt # if true reorient microbe
        turn!(microbe, model)
    end
    nothing
end

# exposed to allow overload and customization
"""
    turnrate(microbe, model)
Evaluate instantaneous turn rate of `microbe`.
"""
turnrate(microbe::AbstractMicrobe, model) = _turnrate(microbe, model)
"""
    affect!(microbe, model)
Can be used to arbitrarily update `microbe` state.
"""
affect!(microbe::AbstractMicrobe, model) = _affect!(microbe, model)
# not exposed and can be used to restore default functionality
_turnrate(microbe::AbstractMicrobe, model) = microbe.turn_rate
_affect!(microbe::AbstractMicrobe, model) = nothing
