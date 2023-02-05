export AbstractXie, Xie, XieNoisy

abstract type AbstractXie{D} <: AbstractMicrobe{D} end

mutable struct Xie{D} <: AbstractXie{D}
    id::Int
    pos::NTuple{D,Float64}
    motility::AbstractMotility
    vel::NTuple{D,Float64}
    speed::Float64
    turn_rate_forward::Float64
    turn_rate_backward::Float64
    rotational_diffusivity::Float64
    radius::Float64
    state::Float64
    state_m::Float64
    state_z::Float64
    adaptation_time_m::Float64
    adaptation_time_z::Float64
    gain_forward::Float64
    gain_backward::Float64
    binding_affinity::Float64
    chemotactic_precision::Float64

    Xie{D}(
        id::Int = rand(1:typemax(Int32)),
        pos::NTuple{D,<:Real} = ntuple(zero, D);
        motility = RunReverseFlick(speed_forward = [46.5]),
        vel::NTuple{D,<:Real} = rand_vel(D),
        speed::Real = rand_speed(motility),
        turn_rate_forward::Real = 2.3, # 1/s
        turn_rate_backward::Real = 1.9, # 1/s
        rotational_diffusivity::Real = 0.26, # rad²/s
        radius::Real = 0.5, # μm
        state::Real = 0.0, # s
        state_m::Real = 0.0, # s
        state_z::Real = 0.0, # s
        adaptation_time_m::Real = 1.29, # s
        adaptation_time_z::Real = 0.28, # s
        gain_forward::Real = 2.7, # 1/s
        gain_backward::Real = 1.6, # 1/s
        binding_affinity::Real = 0.39, # μM
        chemotactic_precision::Real = 0.0,
    ) where {D} = new{D}(
        id, Float64.(pos), motility, Float64.(vel), Float64(speed),
        Float64(turn_rate_forward), Float64(turn_rate_backward),
        Float64(rotational_diffusivity), Float64(radius),
        Float64(state), Float64(state_m), Float64(state_z),
        Float64(adaptation_time_m), Float64(adaptation_time_z),
        Float64(gain_forward), Float64(gain_backward),
        Float64(binding_affinity), Float64(chemotactic_precision)
    )
end # struct


# needs a different show because it does not have a "turn_rate" field
function Base.show(io::IO, ::MIME"text/plain", m::AbstractXie{D}) where D
    println(io, "$(typeof(m)) with $(typeof(m.motility)) motility")
    println(io, "position (μm): $(r2dig.(m.pos)); velocity (μm/s): $(r2dig.(m.vel.*m.speed))")
    println(io, "average unbiased turn rate (Hz): forward $(r2dig(m.turn_rate_forward)), backward $(r2dig(m.turn_rate_backward))")
    s = setdiff(fieldnames(typeof(m)), [:id, :pos, :motility, :vel, :turn_rate_forward, :turn_rate_backward])
    print(io, "other properties: " * join(s, ", "))
end

function turnrate(microbe::AbstractXie, model)
    if microbe.motility isa AbstractMotilityTwoStep
        return turnrate_twostep(microbe, model)
    else
        return turnrate_onestep(microbe, model)
    end
end
function turnrate_twostep(microbe::AbstractXie, model)
    S = microbe.state
    if microbe.motility.state == Forward
        ν₀ = microbe.turn_rate_forward
        β = microbe.gain_forward
    else
        ν₀ = microbe.turn_rate_backward
        β = microbe.gain_backward
    end
    return ν₀*(1 + β*S)
end
function turnrate_onestep(microbe::AbstractXie, model)
    S = microbe.state
    ν₀ = microbe.turn_rate_forward
    β = microbe.gain_forward
    return ν₀*(1 + β*S)
end

function affect!(microbe::AbstractXie, model; ε=1e-16)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    c = model.concentration_field(microbe.pos, model)
    K = microbe.binding_affinity
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = Π * 0.04075 * sqrt(3*c / (5*π*Dc*a*Δt))
    M = rand(Normal(c,σ))
    ϕ = log(1.0 + max(M/K, -1+ε))
    τ_m = microbe.adaptation_time_m
    τ_z = microbe.adaptation_time_z
    a₀ = (τ_m*τ_z)/(τ_m - τ_z)
    m = microbe.state_m
    z = microbe.state_z
    m += (ϕ - m/τ_m)*Δt
    z += (ϕ - z/τ_z)*Δt
    microbe.state_m = m
    microbe.state_z = z
    microbe.state = a₀ * (m/τ_m - z/τ_z)
    return nothing
end