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
mutable struct BrownBerg{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility 
    vel::NTuple{D,Float64}
    speed::Float64
    turn_rate::Float64
    rotational_diffusivity::Float64
    radius::Float64
    state::Float64
    gain::Float64
    receptor_binding_constant::Float64
    memory::Float64

    BrownBerg{D}(
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility::AbstractMotility = RunTumble(),
        vel::NTuple{D,<:Real} = rand_vel(D),
        speed::Real = rand_speed(motility),
        turn_rate::Real = 1/0.67,
        rotational_diffusivity::Real = 0.035,
        radius::Real = 0.5,
        state::Real = 0.0,
        gain::Real = 660.0,
        receptor_binding_constant::Real = 100.0,
        memory::Real = 1.0,
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed), Float64(turn_rate),
        Float64(rotational_diffusivity), Float64(radius), Float64(state),
        Float64(gain), Float64(receptor_binding_constant),
        Float64(memory)
    )
end # struct

function _affect!(microbe::BrownBerg, model)
    Δt = model.timestep
    τₘ = microbe.memory
    β = Δt / τₘ # memory loss factor
    KD = microbe.receptor_binding_constant
    S = microbe.state # weighted dPb/dt at previous step
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    du_dt = dot(microbe.vel, ∇u)*microbe.speed + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    microbe.state = β*M + S*exp(-β) # new weighted dPb/dt
    return nothing
end # function

function _turnrate(microbe::BrownBerg, model)
    ν₀ = microbe.turn_rate # unbiased
    g = microbe.gain
    S = microbe.state
    return ν₀*exp(-g*S) # modulated turn rate
end # function