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
- `motility = RunTumble(speed = [30.0])`
- `turn_rate = 1.49` Hz
- `state = 1`
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
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility = RunTumble(speed = [30.0]),
        vel::NTuple{D,<:Real} = rand_vel(D),
        speed::Real = rand_speed(motility),
        turn_rate::Real = 1/0.67, # 1/s
        markovian_variables::Vector{<:Real} = zeros(3),
        state::Real = 1.0,
        rotational_diffusivity::Real = 0.26, # rad²/s
        gain::Real = 50.0, # 1
        memory::Real = 1.0, # s
        radius::Real = 0.5, # μm
        chemotactic_precision::Real = 0.0,
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed),
        Float64(turn_rate), Float64.(markovian_variables),
        Float64(state), Float64(rotational_diffusivity),
        Float64(gain), Float64(memory), Float64(radius),
        Float64(chemotactic_precision)
    )
end # struct

function _affect!(microbe::Celani, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    c = model.concentration_field(microbe.pos, model)
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = Π * 0.04075 * sqrt(3*c / (5*π*Dc*a*Δt)) # noise (Berg-Purcell)
    M = rand(Normal(c,σ)) # measurement
    λ = 1/microbe.memory
    W = microbe.markovian_variables
    S = microbe.state
    W[1] += (-λ*W[1] + M)*Δt
    W[2] += (-λ*W[2] + W[1])*Δt
    W[3] += (-λ*W[3] + 2*W[2])*Δt
    S = λ^2*(W[2] - λ*W[3]/2)
    return nothing
end # function

function _turnrate(microbe::Celani, model)
    ν₀ = microbe.turn_rate # unbiased
    β = microbe.gain
    S = microbe.state
    return ν₀*(1-β*S) # modulated turn rate
end # function