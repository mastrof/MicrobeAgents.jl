export Celani

"""
    Celani(; gain=50, memory=1, chemotactic_precision=0, field=:chemicalfield)
Chemotaxis behavior using the response kernel from 'Celani and Vergassola (2010) PNAS',
extracted from experiments on E. coli.
Optional sensing noise follows the Berg-Purcell formula, scaled by
`chemotactic_precision` as in 'Brumley et al. (2019) PNAS'; it requires a
microbe `radius > 0`.
Internal markovian variables are initialized at steady state with the local
concentration when the microbe is added.

Parameters:
- `gain = 50`
- `memory = 1` s
- `chemotactic_precision = 0`
- `field = :chemicalfield`: key of the chemical field to sense
Internal state: `state`, `markovian_variables`.
"""
@kwdef mutable struct Celani
    gain::Float64 = 50.0
    memory::Float64 = 1.0
    chemotactic_precision::Float64 = 0.0
    markovian_variables::Vector{Float64} = zeros(3)
    state::Float64 = 0.0
    field::Symbol = :chemicalfield
end

function initialize!(b::Celani, microbe, model)
    check_field(model, b.field)
    check_sensing_radius(b, b.chemotactic_precision, microbe)
    W = b.markovian_variables
    λ = 1 / b.memory
    M = concentration(microbe, model, b.field)
    W[1] = M / λ
    W[2] = W[1] / λ
    W[3] = 2W[2] / λ
    return nothing
end

function affect!(b::Celani, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = diffusivity(microbe, model, b.field)
    c = concentration(microbe, model, b.field)
    a = radius(microbe)
    Π = b.chemotactic_precision
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt)) # noise (Berg-Purcell)
    M = rand(abmrng(model), Normal(c, σ)) # measurement
    λ = 1 / b.memory
    W = b.markovian_variables
    W[1] += (-λ * W[1] + M) * Δt
    W[2] += (-λ * W[2] + W[1]) * Δt
    W[3] += (-λ * W[3] + 2 * W[2]) * Δt
    b.state = λ^2 * (W[2] - λ * W[3] / 2)
    return nothing
end

bias(b::Celani, microbe, model) = 1 - b.gain * b.state
