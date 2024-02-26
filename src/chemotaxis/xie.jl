export Xie

"""
    Xie{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium adapted from 'Xie et al. (2019) Biophys J'.
The model is developed based on experimental measurements of the
chemotactic response function in the marine bacterium V. alginolyticus.
The peculiarity of the model is the presence of distinct parameters for the
forward and backward swimming states.

Sensing noise (not present in the original model) is customarily introduced
through the molecular counting noise formula by Berg and Purcell, and can be
tuned through a `chemotactic_precision` factor inspired by
'Brumley et al. (2019) PNAS' (defaults to 0, i.e. no noise).

Default parameters:
- `motility = RunReverseFlick(speed = [46.5])`
- `turn_rate_forward = 2.3` Hz
- `turn_rate_backward = 1.9` Hz
- `state = 0.0` s
- `state_m = 0.0` s
- `state_z = 0.0` s
- `rotational_diffusivity = 0.26` rad²/s
- `adaptation_time_m = 1.29` s
- `adaptation_time_z = 0.28` s
- `gain_forward = 2.7` 1/s
- `gain_backward = 1.6` 1/s
- `binding_affinity = 0.39` μM
- `chemotactic_precision = 0.0`
- `radius = 0.5` μm
"""
@agent struct Xie{D}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D}
    speed::Float64
    motility = RunReverseFlick(speed = [46.5])
    turn_rate_forward::Float64 = 2.3
    turn_rate_backward::Float64 = 1.9
    rotational_diffusivity::Float64 = 0.26
    radius::Float64 = 0.5
    state::Float64 = 0.0
    state_m::Float64 = 0.0
    state_z::Float64 = 0.0
    adaptation_time_m::Float64 = 1.29
    adaptation_time_z::Float64 = 0.28
    gain_forward::Float64 = 2.7
    gain_backward::Float64 = 1.6
    binding_affinity::Float64 = 0.39
    chemotactic_precision::Float64 = 0.0
end

# needs a different show because it does not have a "turn_rate" field
function Base.show(io::IO, ::MIME"text/plain", m::Xie{D}) where {D}
    println(io, "$(typeof(m)) with $(typeof(m.motility)) motility")
    println(io, "position (μm): $(r2dig.(m.pos)); velocity (μm/s): $(r2dig.(m.vel.*m.speed))")
    println(io, "average unbiased turn rate (Hz): forward $(r2dig(m.turn_rate_forward)), backward $(r2dig(m.turn_rate_backward))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate_forward, :turn_rate_backward])
    print(io, "other properties: " * join(s, ", "))
end

function chemotaxis!(microbe::Xie, model; ε=1e-16)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    c = model.concentration_field(microbe.pos, model)
    K = microbe.binding_affinity
    a = microbe.radius
    Π = microbe.chemotactic_precision
    # noisy concentration measurement with Berg-Purcell formula
    σ = CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt))
    M = rand(Normal(c, σ))
    ϕ = log(1.0 + max(M / K, -1 + ε))
    τ_m = microbe.adaptation_time_m
    τ_z = microbe.adaptation_time_z
    a₀ = (τ_m * τ_z) / (τ_m - τ_z)
    m = microbe.state_m
    z = microbe.state_z
    m += (ϕ - m / τ_m) * Δt
    z += (ϕ - z / τ_z) * Δt
    microbe.state_m = m
    microbe.state_z = z
    microbe.state = a₀ * (m / τ_m - z / τ_z)
    return nothing
end

function tumblebias(microbe::Xie)
    S = state(microbe)
    if state(motilepattern(microbe)) == Forward
        β = microbe.gain_forward
    else
        β = microbe.gain_backward
    end
    return (1 + β*S)
end

function turnrate(microbe::Xie)
    if state(motilepattern(microbe)) == Forward
        return microbe.turn_rate_forward
    else
        return microbe.turn_rate_backward
    end
end
