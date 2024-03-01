export BrownBerg

"""
    BrownBerg{D} <: AbstractMicrobe{D}
Model of chemotactic E.coli from 'Brown and Berg (1974) PNAS'

Default parameters:
- `motility = RunTumble()`
- `turn_rate = 1.49` Hz frequency of reorientations
- `rotational_diffusivity = 0.035` rad²/s coefficient of brownian rotational diffusion
- `radius = 0.5` μm equivalent spherical radius of the microbe
- `state = 0.0` corresponds to 'weighted dPb/dt' in the paper
- `gain = 660` s
- `receptor_binding_constant = 100` μM
- `memory = 1` s
"""
@agent struct BrownBerg{D}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D}
    speed::Float64
    motility = RunTumble()
    turn_rate::Float64 = 1 / 0.67
    rotational_diffusivity::Float64 = 0.035
    radius::Float64 = 0.5
    state::Float64 = 0.0
    gain::Float64 = 660.0
    receptor_binding_constant::Float64 = 100.0
    memory::Float64 = 1.0
end

function chemotaxis!(microbe::BrownBerg, model)
    Δt = model.timestep
    τₘ = microbe.memory
    β = exp(-Δt / τₘ) # memory loss factor
    KD = microbe.receptor_binding_constant
    S = state(microbe) # weighted dPb/dt at previous step
    pos = position(microbe)
    vel = velocity(microbe)
    u = concentration(pos, model)
    ∇u = gradient(pos, model)
    ∂ₜu = time_derivative(pos, model)
    du_dt = dot(vel, ∇u) + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    microbe.state = (1 - β) * M + S * β # new weighted dPb/dt
    return nothing
end # function

function tumblebias(microbe::BrownBerg)
    g = microbe.gain
    S = state(microbe)
    return exp(-g*S)
end
