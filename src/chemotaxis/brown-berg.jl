export BrownBerg

"""
    BrownBerg(; gain=660, receptor_binding_constant=100, memory=1, field=:chemicalfield)
Chemotaxis behavior from 'Brown and Berg (1974) PNAS'.

Parameters:
- `gain = 660` s
- `receptor_binding_constant = 100` μM
- `memory = 1` s
- `field = :chemicalfield`: key of the chemical field to sense
Internal state: `state`, the weighted dPb/dt of the paper.
"""
@kwdef mutable struct BrownBerg
    gain::Float64 = 660.0
    receptor_binding_constant::Float64 = 100.0
    memory::Float64 = 1.0
    state::Float64 = 0.0
    field::Symbol = :chemicalfield
end

initialize!(b::BrownBerg, microbe, model) = check_field(model, b.field)

function affect!(b::BrownBerg, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    τₘ = b.memory
    β = exp(-Δt / τₘ) # memory loss factor
    KD = b.receptor_binding_constant
    S = b.state # weighted dPb/dt at previous step
    vel = velocity(microbe)
    u = concentration(microbe, model, b.field)
    ∇u = gradient(microbe, model, b.field)
    ∂ₜu = time_derivative(microbe, model, b.field)
    du_dt = dot(vel, ∇u) + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    b.state = (1 - β) * M + S * β # new weighted dPb/dt
    return nothing
end

bias(b::BrownBerg, microbe, model) = exp(-b.gain * b.state)
