export Brumley, chemotaxis!

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
- `memory = 1.3` s → 'τₘ'
- `gain_receptor = 50.0` μM⁻¹ → 'κ'
- `gain = 50.0` → 'Γ'
- `chemotactic_precision = 6.0` → 'Π'
- `radius = 0.5` μm → 'a'
"""
@agent Brumley{D} ContinuousAgent{D} where D AbstractMicrobe{D} begin
    speed::Float64
    motility::AbstractMotility = RunReverseFlick(speed_forward = [46.5]),
    turn_rate::Float64 = 1 / 0.45
    rotational_diffusivity::Float64 = 0.035
    radius::Float64 = 0.5
    state::Float64 = 0.0
    memory::Float64 = 1.3
    gain_receptor::Float64 = 50.0
    gain::Float64 = 50.0
    chemotactic_precision::Float64 = 6.0
end

function chemotaxis!(microbe::Brumley, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    τₘ = microbe.memory
    α = exp(-Δt / τₘ) # memory persistence factor
    a = microbe.radius
    Π = microbe.chemotactic_precision
    κ = microbe.gain_receptor
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    # gradient measurement
    μ = dot(microbe.vel, ∇u) * microbe.speed + ∂ₜu # mean
    σ = CONV_NOISE * Π * sqrt(3 * u / (π * a * Dc * Δt^3)) # noise
    M = rand(Normal(μ, σ)) # measurement
    # update internal state
    S = microbe.state
    microbe.state = α * S + (1 - α) * κ * τₘ * M
    return nothing
end # function

function affect!(microbe::Brumley, model)
    chemotaxis!(microbe, model)
end

function turnrate(microbe::Brumley, model)
    ν₀ = microbe.turn_rate # unbiased
    Γ = microbe.gain
    S = microbe.state
    return (1 + exp(-Γ * S)) * ν₀ / 2 # modulated turn rate
end # function
