export SonMenolascina

@agent struct SonMenolascina{D,N}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N}
    rotational_diffusivity::Float64 = 0.035
    radius::Float64 = 0.5
    state::Float64 = 0.0
    gain::Float64 = 660.0
    receptor_binding_constant::Float64 = 100.0
    memory::Float64 = 1.0
    #= values from paper -- not super good and units are wrong
    eta::Float64 = -0.3942 # s/μm
    ζ::Float64 = -0.2019 # s/μm
    vT::Float64 = 18.88 # μm/s
    θ::Float64 = 0.8452 # s/μm
    =#
    eta::Float64 = -0.55 # s
    ζ::Float64 = -0.35 # s/μm
    vT::Float64 = 18.88 # μm/s
    θ::Float64 = 1.0 # s
    cT::Float64 = 0.05 # μM
    chemokinesis_on::Bool = false # turned on when local c > cT
    chemokinetic_factor::Float64 = 1.3
end

function speed(microbe::SonMenolascina)
    if microbe.chemokinesis_on
        return microbe.speed * microbe.chemokinetic_factor
    else
        return microbe.speed
    end
end

function chemokinesis!(microbe::SonMenolascina, c)
    microbe.chemokinesis_on = (c >= microbe.cT)
end

function chemotaxis!(microbe::SonMenolascina, model)
    Δt = abmtimestep(model)
    τₘ = microbe.memory
    β = exp(-Δt / τₘ) # memory loss factor
    KD = microbe.receptor_binding_constant
    S = state(microbe) # weighted dPb/dt at previous step
    pos = position(microbe)
    u = concentration(pos, model)
    chemokinesis!(microbe, u) # modify speed based on local concentration
    ∇u = gradient(pos, model)
    ∂ₜu = time_derivative(pos, model)
    vel = velocity(microbe)
    du_dt = dot(vel, ∇u) + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    microbe.state = (1 - β) * M + S * β # new weighted dPb/dt
    return nothing
end # function

# speed-dependent flick probability if 4-state motility is used
function update_motilestate!(microbe::SonMenolascina{D,4}, model::AgentBasedModel) where D
    motility = motilepattern(microbe)
    i = state(motility)
    w = transition_weights(motility, i)
    if i == 3 # backward run
        p_flick = flick_probability(microbe, model)
        w[2] = 1-p_flick
        w[4] = p_flick
    end
    j = sample(abmrng(model), eachindex(w), w)
    update_motilestate!(motility, j)
end
# turn hard-coded values from paper into parameters?
function flick_probability(microbe::SonMenolascina, model)
    0.055 + 0.72 / (1 + exp(-0.25*(speed(microbe)-36.0)))
end

function switching_probability(microbe::SonMenolascina, model)
    dt = abmtimestep(model)
    M = motilestate(microbe)
    τ = duration(M)
    β = biased(M) ? bias(microbe)*speed_dependent_bias(microbe) : 1.0
    return β * dt / τ
end

function bias(microbe::SonMenolascina)
    g = microbe.gain
    S = state(microbe)
    return exp(-g*S)
end

function speed_dependent_bias(microbe::SonMenolascina)
    n = microbe.eta
    θ = microbe.θ
    ζ = microbe.ζ
    vT = microbe.vT
    v = speed(microbe)
    f = speed_dependent_turnrate(v, n, ζ, θ, vT)
    f0 = speed_dependent_turnrate(zero(v), n, ζ, θ, vT)
    f / f0
end

function speed_dependent_turnrate(v, n, ζ, θ, vT)
    1 / (n / (1 + exp(ζ*(v-vT))) + θ)
end
