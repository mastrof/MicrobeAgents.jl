export Xie

"""
    Xie(; adaptation_time_m=1.29, adaptation_time_z=0.28,
        gain_forward=2.7, gain_backward=1.6, binding_affinity=0.39,
        chemotactic_precision=0)
Chemotaxis behavior adapted from 'Xie et al. (2019) Biophys J', based on the
response function of V. alginolyticus. With a 4-state motility
(`RunReverseFlick`), `gain_backward` applies in the backward run (state 3)
and `gain_forward` otherwise.
Optional Berg-Purcell sensing noise (`chemotactic_precision`) requires a
microbe `radius > 0`.

Parameters:
- `adaptation_time_m = 1.29` s
- `adaptation_time_z = 0.28` s
- `gain_forward = 2.7` 1/s
- `gain_backward = 1.6` 1/s
- `binding_affinity = 0.39` μM
- `chemotactic_precision = 0`
Internal state: `state`, `state_m`, `state_z`.
"""
@kwdef mutable struct Xie
    adaptation_time_m::Float64 = 1.29
    adaptation_time_z::Float64 = 0.28
    gain_forward::Float64 = 2.7
    gain_backward::Float64 = 1.6
    binding_affinity::Float64 = 0.39
    chemotactic_precision::Float64 = 0.0
    state::Float64 = 0.0
    state_m::Float64 = 0.0
    state_z::Float64 = 0.0
end

initialize!(b::Xie, microbe, model) =
    check_sensing_radius(b, b.chemotactic_precision, microbe)

function affect!(b::Xie, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = diffusivity(microbe, model)
    c = concentration(microbe, model)
    K = b.binding_affinity
    a = radius(microbe)
    Π = b.chemotactic_precision
    # noisy concentration measurement with Berg-Purcell formula
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt))
    M = max(rand(abmrng(model), Normal(c, σ)), zero(c))
    ϕ = log(1.0 + M / K)
    τ_m = b.adaptation_time_m
    τ_z = b.adaptation_time_z
    a₀ = (τ_m * τ_z) / (τ_m - τ_z)
    m = b.state_m
    z = b.state_z
    m += (ϕ - m / τ_m) * Δt
    z += (ϕ - z / τ_z) * Δt
    b.state_m = m
    b.state_z = z
    b.state = a₀ * (m / τ_m - z / τ_z)
    return nothing
end

function bias(b::Xie, microbe::AbstractMicrobe{D,N}, model) where {D,N}
    backward = N == 4 && state(motilepattern(microbe)) == 3
    β = backward ? b.gain_backward : b.gain_forward
    return 1 + β * b.state
end
