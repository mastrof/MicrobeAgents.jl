export Celani, chemotaxis!

"""
    Celani{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium using the response kernel from
'Celani and Vergassola (2010) PNAS', extracted from experiments on E. coli.

Sensing noise (not present in the original model) is customarily introduced
through the molecular counting noise formula by Berg and Purcell, and can be
tuned through a `chemotactic_precision` factor inspired by
'Brumley et al. (2019) PNAS' (defaults to 0, i.e. no noise).

Default parameters:
- `motility = RunTumble(speed = [30.0])`
- `turn_rate = 1.49` Hz
- `state = 0`
- `rotational_diffusivity = 0.26` rad²/s
- `gain = 50.0`
- `memory = 1` s
- `radius = 0.5` μm
"""
mutable struct Celani{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility
    vel::NTuple{D,Float64}
    speed::Float64
    turn_rate::Float64
    markovian_variables::Vector{Float64}
    state::Float64
    rotational_diffusivity::Float64
    gain::Float64
    memory::Float64
    radius::Float64
    chemotactic_precision::Float64

    Celani{D}(
        id::Int=rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real}=ntuple(zero, D);
        motility=RunTumble(speed=[30.0]),
        vel::NTuple{D,<:Real}=rand_vel(D),
        speed::Real=rand_speed(motility),
        turn_rate::Real=1 / 0.67, # 1/s
        markovian_variables::Vector{<:Real}=zeros(3),
        state::Real=0.0,
        rotational_diffusivity::Real=0.26, # rad²/s
        gain::Real=50.0, # 1
        memory::Real=1.0, # s
        radius::Real=0.5, # μm
        chemotactic_precision::Real=0.0
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed),
        Float64(turn_rate), Float64.(markovian_variables),
        Float64(state), Float64(rotational_diffusivity),
        Float64(gain), Float64(memory), Float64(radius),
        Float64(chemotactic_precision)
    )
end # struct

function chemotaxis!(microbe::Celani, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    c = model.concentration_field(microbe.pos, model)
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt)) # noise (Berg-Purcell)
    M = rand(Normal(c, σ)) # measurement
    λ = 1 / microbe.memory
    W = microbe.markovian_variables
    W[1] += (-λ * W[1] + M) * Δt
    W[2] += (-λ * W[2] + W[1]) * Δt
    W[3] += (-λ * W[3] + 2 * W[2]) * Δt
    microbe.state = λ^2 * (W[2] - λ * W[3] / 2)
    return nothing
end # function

function affect!(microbe::Celani, model)
    chemotaxis!(microbe, model)
end

function turnrate(microbe::Celani, model)
    ν₀ = microbe.turn_rate # unbiased
    β = microbe.gain
    S = microbe.state
    return ν₀ * (1 - β * S) # modulated turn rate
end # function

# Celani requires a custom add_agent! method
# to initialize the markovian variables at steady state
# depending the concentration field at the agent position
function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{Celani{D}},
    model::AgentBasedModel,
    properties...;
    vel=nothing,
    speed=nothing,
    kwargs...
) where {D}
    id = nextid(model)
    microbe = A(id, pos, properties...; vel=ntuple(zero, D), speed=0, kwargs...)
    microbe.vel = isnothing(vel) ? rand_vel(abmrng(model), D) : vel
    microbe.speed = isnothing(speed) ? rand_speed(abmrng(model), microbe.motility) : speed
    initialize_markovian_variables!(microbe, model)
    add_agent_pos!(microbe, model)
end

"""
    initialize_markovian_variables!(microbe::Celani, model)
Initialize the internal markovian variables of `microbe` to be at steady state
with respect to the `concentration_field` defined by `model`.
"""
function initialize_markovian_variables!(microbe::Celani, model)
    W = microbe.markovian_variables
    λ = 1 / microbe.memory
    M = model.concentration_field(microbe.pos, model)
    W[1] = M / λ
    W[2] = W[1] / λ
    W[3] = 2W[2] / λ
    nothing
end
