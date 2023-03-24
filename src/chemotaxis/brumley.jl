export Brumley

"""
    Brumley{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium from 'Brumley et al. (2019) PNAS'.
The model is optimized for simulation of marine bacteria and accounts
for the presence of (gaussian) sensing noise in the chemotactic pathway.

Default parameters:
- `motility = RunReverseFlick(speed_forward = [46.5])`
- `turn_rate = 2.22` Hz → '1/τ₀'
- `state = 0.0` → 'S'
- `rotational_diffusivity = 0.035` rad²/s
- `adaptation_time = 1.3` s → 'τₘ'
- `receptor_gain = 50.0` μM⁻¹ → 'κ'
- `motor_gain = 50.0` → 'Γ'
- `chemotactic_precision = 6.0` → 'Π'
- `radius = 0.5` μm → 'a'
"""
mutable struct Brumley{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility
    vel::NTuple{D,Float64}
    speed::Float64
    turn_rate::Float64
    rotational_diffusivity::Float64
    radius::Float64
    state::Float64
    adaptation_time::Float64
    receptor_gain::Float64
    motor_gain::Float64
    chemotactic_precision::Float64

    Brumley{D}(
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility = RunReverseFlick(speed_forward = [46.5]),
        vel::NTuple{D,<:Real} = rand_vel(D),
        speed::Real = rand_speed(motility),
        turn_rate::Real = 1/0.45,
        rotational_diffusivity::Real = 0.035,
        radius::Real = 0.5,
        state::Real = 0.0,
        adaptation_time::Real = 1.3,
        receptor_gain::Real = 50.0,
        motor_gain::Real = 50.0,
        chemotactic_precision::Real = 6.0,
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed), Float64(turn_rate),
        Float64(rotational_diffusivity), Float64(radius), Float64(state), 
        Float64(adaptation_time), Float64(receptor_gain),
        Float64(motor_gain), Float64(chemotactic_precision)
    )

end # struct

function _affect!(microbe::Brumley, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    τₘ = microbe.adaptation_time
    α = exp(-Δt/τₘ) # memory persistence factor
    a = microbe.radius
    Π = microbe.chemotactic_precision
    κ = microbe.receptor_gain
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    # gradient measurement
    μ = dot(microbe.vel, ∇u)*microbe.speed + ∂ₜu # mean
    σ = Π * 0.04075 * sqrt(3*u / (π*a*Dc*Δt^3)) # noise
    M = rand(Normal(μ,σ)) # measurement
    # update internal state
    S = microbe.state
    microbe.state = α*S + (1-α)*κ*τₘ*M
    return nothing
end # function

function turnrate(microbe::Brumley, model)
    ν₀ = microbe.turn_rate # unbiased
    Γ = microbe.motor_gain
    S = microbe.state
    return (1 + exp(-Γ*S)) * ν₀/2 # modulated turn rate
end # function
