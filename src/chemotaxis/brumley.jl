export Brumley

"""
    Brumley{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium from 'Brumley et al. (2019) PNAS'.
The model is optimized for simulation of marine bacteria and accounts
for the presence of (gaussian) sensing noise in the chemotactic pathway.

Default parameters:
- `motility = RunReverseFlick(0.45, [46.5], 0.45, [46.5])`
- `state = 0.0` → 'S'
- `rotational_diffusivity = 0.035` rad²/s
- `memory = 1.3` s → 'τₘ'
- `gain_receptor = 50.0` μM⁻¹ → 'κ'
- `gain = 50.0` → 'Γ'
- `chemotactic_precision = 6.0` → 'Π'
- `radius = 0.5` μm → 'a'
"""
@agent struct Brumley{D,N}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N} = RunReverseFlick(;
        run_speed_forward=[46.5], run_duration_forward=0.45,
        run_speed_backward=[46.5], run_duration_backward=0.45,
    )
    rotational_diffusivity::Float64 = 0.035
    radius::Float64 = 0.5
    state::Float64 = 0.0
    memory::Float64 = 1.3
    gain_receptor::Float64 = 50.0
    gain::Float64 = 50.0
    chemotactic_precision::Float64 = 6.0
end

function chemotaxis!(microbe::Brumley, model)
    Δt = abmtimestep(model)
    Dc = chemoattractant_diffusivity(model)
    τₘ = microbe.memory
    α = exp(-Δt / τₘ) # memory persistence factor
    a = microbe.radius
    Π = microbe.chemotactic_precision
    κ = microbe.gain_receptor
    pos = position(microbe)
    vel = velocity(microbe)
    u = concentration(pos, model)
    ∇u = gradient(pos, model)
    ∂ₜu = time_derivative(pos, model)
    # gradient measurement
    μ = dot(vel, ∇u) + ∂ₜu # mean
    σ = CONV_NOISE * Π * sqrt(3 * u / (π * a * Dc * Δt^3)) # noise
    M = rand(abmrng(model), Normal(μ, σ)) # measurement
    # update internal state
    S = state(microbe)
    microbe.state = α * S + (1 - α) * κ * τₘ * M
    return nothing
end # function

function bias(microbe::Brumley)
    Γ = microbe.gain
    S = state(microbe)
    return (1 + exp(-Γ*S))/2
end
