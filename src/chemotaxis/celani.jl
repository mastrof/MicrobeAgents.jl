export Celani

"""
    Celani{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium using the response kernel from
'Celani and Vergassola (2010) PNAS', extracted from experiments on E. coli.

Sensing noise (not present in the original model) is customarily introduced
through the molecular counting noise formula by Berg and Purcell, and can be
tuned through a `chemotactic_precision` factor inspired by
'Brumley et al. (2019) PNAS' (defaults to 0, i.e. no noise).

Default parameters:
- `motility = RunTumble(0.67, [30.0], 0.1)
- `state = 0`
- `rotational_diffusivity = 0.26` rad²/s
- `gain = 50.0`
- `memory = 1` s
- `radius = 0.5` μm
"""
@agent struct Celani{D,N}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N} = RunTumble(0.67, [30.0], 0.1)
    rotational_diffusivity::Float64 = 0.26
    radius::Float64 = 0.5
    state::Float64 = 0.0
    markovian_variables::Vector{Float64} = zeros(3)
    gain::Float64 = 50.0
    memory::Float64 = 1.0
    chemotactic_precision::Float64 = 0.0
end

function chemotaxis!(microbe::Celani, model)
    Δt = abmtimestep(model)
    Dc = chemoattractant_diffusivity(model)
    c = concentration(position(microbe), model)
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt)) # noise (Berg-Purcell)
    M = rand(abmrng(model), Normal(c, σ)) # measurement
    λ = 1 / microbe.memory
    W = microbe.markovian_variables
    W[1] += (-λ * W[1] + M) * Δt
    W[2] += (-λ * W[2] + W[1]) * Δt
    W[3] += (-λ * W[3] + 2 * W[2]) * Δt
    microbe.state = λ^2 * (W[2] - λ * W[3] / 2)
    return nothing
end # function

function bias(microbe::Celani)
    β = microbe.gain
    S = state(microbe)
    return (1 - β*S)
end

# Celani requires a custom add_agent! method
# to initialize the markovian variables at steady state
# depending on the concentration field at the agent position
function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{<:Celani{D}},
    model::AgentBasedModel,
    properties...;
    vel=nothing,
    speed=nothing,
    kwproperties...
) where {D}
    @assert haskey(kwproperties, :motility) "Missing required keyword argument `motility`"
    N = get_motility_N(kwproperties[:motility])
    add_agent!(pos, A{N}, model, properties...; vel, speed, kwproperties...)
end

function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{Celani{D,N}},
    model::AgentBasedModel,
    properties...;
    vel=nothing,
    speed=nothing,
    kwproperties...
) where {D,N}
    @assert haskey(kwproperties, :motility) "Missing required keyword argument `motility`"
    id = Agents.nextid(model)
    if !isempty(properties)
        microbe = A(id, pos, properties...)
    else
        microbe = A(; id, pos, vel = zero(SVector{D}), speed = 0.0, kwproperties...)
        microbe.vel = isnothing(vel) ? random_velocity(model) : vel
        microbe.speed = isnothing(speed) ? random_speed(microbe, model) : speed
        initialize_markovian_variables!(microbe, model)
    end
    Agents.add_agent_own_pos!(microbe, model)
end

"""
    initialize_markovian_variables!(microbe::Celani, model)
Initialize the internal markovian variables of `microbe` to be at steady state
with respect to the `concentration_field` defined by `model`.
"""
function initialize_markovian_variables!(microbe::Celani, model)
    W = microbe.markovian_variables
    λ = 1 / microbe.memory
    M = concentration(position(microbe), model)
    W[1] = M / λ
    W[2] = W[1] / λ
    W[3] = 2W[2] / λ
    nothing
end
