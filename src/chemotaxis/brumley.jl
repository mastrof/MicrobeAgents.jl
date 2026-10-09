export Brumley

"""
    Brumley(; memory=1.3, gain_receptor=50, gain=50, chemotactic_precision=6)
Chemotaxis behavior from 'Brumley et al. (2019) PNAS', with gaussian sensing
noise in the gradient measurement. Requires a microbe `radius > 0` when
`chemotactic_precision > 0`.

Parameters:
- `memory = 1.3` s → 'τₘ'
- `gain_receptor = 50` μM⁻¹ → 'κ'
- `gain = 50` → 'Γ'
- `chemotactic_precision = 6` → 'Π'
Internal state: `state` → 'S'.
"""
@kwdef mutable struct Brumley
    memory::Float64 = 1.3
    gain_receptor::Float64 = 50.0
    gain::Float64 = 50.0
    chemotactic_precision::Float64 = 6.0
    state::Float64 = 0.0
end

initialize!(b::Brumley, microbe, model) =
    check_sensing_radius(b, b.chemotactic_precision, microbe)

function affect!(b::Brumley, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = diffusivity(microbe, model)
    τₘ = b.memory
    α = exp(-Δt / τₘ) # memory persistence factor
    a = radius(microbe)
    Π = b.chemotactic_precision
    κ = b.gain_receptor
    vel = velocity(microbe)
    u = concentration(microbe, model)
    ∇u = gradient(microbe, model)
    ∂ₜu = time_derivative(microbe, model)
    # gradient measurement
    μ = dot(vel, ∇u) + ∂ₜu # mean
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * u / (π * a * Dc * Δt^3)) # noise
    M = rand(abmrng(model), Normal(μ, σ)) # measurement
    # update internal state
    S = b.state
    b.state = α * S + (1 - α) * κ * τₘ * M
    return nothing
end

bias(b::Brumley, microbe, model) = (1 + exp(-b.gain * b.state)) / 2
