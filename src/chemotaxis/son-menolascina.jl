export SonMenolascina, Chemokinesis, SpeedDependentTurnRate, SpeedDependentFlick

"""
    Chemokinesis(; threshold=0.05, factor=1.3, field=:chemicalfield)
Multiply the microbe speed by `factor` while the local concentration of the
chemical field `field` is at least `threshold` (μM).
"""
@kwdef mutable struct Chemokinesis
    threshold::Float64 = 0.05
    factor::Float64 = 1.3
    on::Bool = false
    field::Symbol = :chemicalfield
end

initialize!(b::Chemokinesis, microbe, model) = check_field(model, b.field)

function affect!(b::Chemokinesis, microbe::AbstractMicrobe, model)
    b.on = concentration(microbe, model, b.field) >= b.threshold
    return nothing
end
speed_factor(b::Chemokinesis, microbe) = b.on ? b.factor : 1.0

"""
    SpeedDependentTurnRate(; eta=-0.55, ζ=-0.35, θ=1.0, vT=18.88)
Bias the switching rate by `f(v)/f(0)` with
`f(v) = 1 / (eta / (1 + exp(ζ*(v - vT))) + θ)` and `v = speed(microbe)`
(from 'Son, Menolascina and Stocker (2016) PNAS').
"""
@kwdef struct SpeedDependentTurnRate
    eta::Float64 = -0.55 # s
    ζ::Float64 = -0.35 # s/μm
    θ::Float64 = 1.0 # s
    vT::Float64 = 18.88 # μm/s
end

_turnrate(v, b::SpeedDependentTurnRate) = 1 / (b.eta / (1 + exp(b.ζ * (v - b.vT))) + b.θ)
function bias(b::SpeedDependentTurnRate, microbe, model)
    v = speed(microbe)
    return _turnrate(v, b) / _turnrate(zero(v), b)
end

"""
    SpeedDependentFlick(; p0=0.055, amplitude=0.72, steepness=0.25, v_half=36.0)
With a 4-state motility (`RunReverseFlick`), make the probability of flicking
after a backward run depend on speed:
`p = p0 + amplitude / (1 + exp(-steepness*(speed - v_half)))`;
otherwise the run ends in a reversal.
"""
@kwdef struct SpeedDependentFlick
    p0::Float64 = 0.055
    amplitude::Float64 = 0.72
    steepness::Float64 = 0.25 # s/μm
    v_half::Float64 = 36.0 # μm/s
end

flick_probability(b::SpeedDependentFlick, microbe) =
    b.p0 + b.amplitude / (1 + exp(-b.steepness * (speed(microbe) - b.v_half)))

function transition_weights!(w, b::SpeedDependentFlick, microbe::AbstractMicrobe{D,N}, model) where {D,N}
    if N == 4 && state(motilepattern(microbe)) == 3 # backward run
        p = flick_probability(b, microbe)
        w[2] = 1 - p
        w[4] = p
    end
    return nothing
end

"""
    SonMenolascina(; gain=660, receptor_binding_constant=100, memory=1,
        eta=-0.55, ζ=-0.35, θ=1.0, vT=18.88, threshold=0.05, factor=1.3,
        field=:chemicalfield)
Behaviors of the chemotaxis model from 'Son, Menolascina and Stocker (2016) PNAS':
returns the `NamedTuple`
`(chemokinesis = Chemokinesis(...), chemotaxis = BrownBerg(...),
turnrate = SpeedDependentTurnRate(...), flick = SpeedDependentFlick())`,
to be passed as `behaviors` (use with a `RunReverseFlick` motility).
`field` is the key of the chemical field sensed by chemokinesis and chemotaxis.
"""
function SonMenolascina(;
    gain = 660.0, receptor_binding_constant = 100.0, memory = 1.0,
    eta = -0.55, ζ = -0.35, θ = 1.0, vT = 18.88,
    threshold = 0.05, factor = 1.3, field = :chemicalfield,
)
    (
        chemokinesis = Chemokinesis(; threshold, factor, field),
        chemotaxis = BrownBerg(; gain, receptor_binding_constant, memory, field),
        turnrate = SpeedDependentTurnRate(; eta, ζ, θ, vT),
        flick = SpeedDependentFlick(),
    )
end
