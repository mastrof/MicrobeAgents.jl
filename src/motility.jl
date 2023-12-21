export
    AbstractMotility, MotilityOneStep, MotilityTwoStep,
    RunTumble, RunReverse, RunReverseFlick,
    MotileState, TwoState, Forward, Backward, switch!,
    state, speed, polar, azimuthal

export Arccos # from Agents

"""
    AbstractMotility
General abstract interface for motility patterns.
"""
abstract type AbstractMotility end

"""
    MotilityOneStep
Type for one-step motility patterns (e.g. `RunTumble`).

A `MotilityOneStep` has the fields
- `speed`: distribution of microbe speed, new values extracted after each turn
- `polar`: distribution of polar angles
- `azimuthal`: distribution azimuthal angles
For 2-dimensional microbe types, only `polar` defines reorientations and `azimuthal` is ignored.
"""
struct MotilityOneStep <: AbstractMotility
    speed
    polar
    azimuthal
end

"""
    RunTumble(; speed=(30.0,), polar=Uniform(-π,π), azimuthal=Arccos())
Constructor for a `MotilityOneStep` with default values associated to run-and-tumble motion.
"""
RunTumble(; speed=(30.0,), polar=Uniform(-π,π), azimuthal=Arccos()) =
    MotilityOneStep(speed, polar, azimuthal)


"""
    MotilityTwoStep
Type for two-step motility patterns (e.g. `RunReverse`, `RunReverseFlick`)
In two-step motility patterns, the two "steps" can have different properties.

A `MotilityTwoStep` has the fields
- `speed`: distribution of microbe speed, new values extracted after each turn
- `polar`: distribution of in-plane reorientations for motile state fw
- `azimuthal`: distribution of out-of-plane reorientations for motile state fw
- `speed_backward`: distribution of microbe speed, new values extracted after each turn
- `polar_backward`: distribution of in-plane reorientations for motile state bw
- `azimuthal_backward`: distribution of out-of-plane reorientations for motile state bw
- `motile_state`: defines current motile state (e.g. `Forward` or `Backward` for a `TwoState`)
For 2-dimensional microbe types, only polar distributions define reorientations
while azimuthal ones are ignored.
"""
struct MotilityTwoStep <: AbstractMotility
    speed
    polar
    azimuthal
    speed_backward
    polar_backward
    azimuthal_backward
    motile_state
end

"""
    RunReverse(; speed=(30.0,), polar=(π,), azimuthal=Arccos(), speed_backward=speed, polar_backward=polar, azimuthal_backward=azimuthal)
Constructor for a `MotilityTwoStep` with default values associated to run-reverse motion.
"""
RunReverse(;
    speed = (30.0,),
    polar = (π,),
    azimuthal = Arccos(),
    speed_backward = speed,
    polar_backward = polar,
    azimuthal_backward = azimuthal,
    motile_state = MotileState()
) = MotilityTwoStep(speed, polar, azimuthal,
                    speed_backward, polar_backward, azimuthal_backward,
                    motile_state)

"""
    RunReverseFlick(; speed=(30.0,), polar=(π,), azimuthal=Arccos(), speed_backward=speed, polar_backward=(-π/2,π/2), azimuthal_backward=azimuthal)
Constructor for a `MotilityTwoStep` with default values associated to run-reverse-flick motion.
"""
RunReverseFlick(;
    speed = (30.0,),
    polar = (π,),
    azimuthal = Arccos(),
    speed_backward = speed,
    polar_backward = (-π/2, π/2),
    azimuthal_backward = azimuthal,
    motile_state = MotileState()
) = MotilityTwoStep(speed, polar, azimuthal,
                    speed_backward, polar_backward, azimuthal_backward,
                    motile_state)



# just a wrapper to allow state to be mutable
mutable struct MotileState
    state
end
MotileState() = MotileState(TwoState())
"""
    TwoState <: Enum{Bool}
Represents the state of a two-step motile pattern, can take values
`Forward` or `Backward`.
"""
@enum TwoState::Bool Forward Backward
Base.show(io::IO, ::MIME"text/plain", x::TwoState) = 
    x == Forward ? print(io, "Forward") : print(io, "Backward")
# initialize to Forward if not specified
TwoState() = Forward
# overload getproperty and setproperty! for convenient access to state
function Base.getproperty(obj::MotilityTwoStep, sym::Symbol)
    if sym === :state
        return obj.motile_state.state
    else
        return getfield(obj, sym)
    end
end
function Base.setproperty!(obj::MotilityTwoStep, sym::Symbol, x)
    if sym === :state
        return setfield!(obj.motile_state, :state, x)
    else
        return setfield!(obj, sym, x)
    end
end
# define rules for switching motile state
switch!(::MotilityOneStep) = nothing
""" 
    switch!(m::MotilityTwoStep)
Switch the state of a two-step motility pattern (`m.state`)
from `Forward` to `Backward` and viceversa.
"""
switch!(m::MotilityTwoStep) = (m.state = ~m.state; nothing)
Base.:~(x::TwoState) = TwoState(~Bool(x))

# convenient wrapper to sample new motile state
Base.rand(rng::AbstractRNG, motility::AbstractMotility) = (
    rand(rng, speed(motility)),
    rand(rng, polar(motility)),
    rand(rng, azimuthal(motility))
)

state(m::MotilityOneStep) = Forward
state(m::MotilityTwoStep) = m.state
speed(m::MotilityOneStep) = m.speed
polar(m::MotilityOneStep) = m.polar
azimuthal(m::MotilityOneStep) = m.azimuthal
speed(m::MotilityTwoStep) = state(m) == Forward ? m.speed : m.speed_backward
polar(m::MotilityTwoStep) = state(m) == Forward ? m.polar : m.polar_backward
azimuthal(m::MotilityTwoStep) = state(m) == Forward ? m.azimuthal : m.azimuthal_backward
