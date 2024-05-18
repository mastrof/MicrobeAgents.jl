export MotileState, RunState, TurnState, Run, Tumble, Reverse, Flick, Stop
export Motility, RunTumble, RunReverse, RunReverseFlick, RunStop
export update_motilestate!, motilestate, state, states, transition_weights, duration
export Arccos # from Agents

@compact_structs MotileState begin
    @kwdef struct RunState
        duration
        speed
        polar = nothing
        azimuthal = nothing
    end
    @kwdef struct TurnState
        duration
        speed = [zero(duration)]
        polar
        azimuthal
    end
end

function biased(s::MotileState)
    if kindof(s) === :RunState
        return true
    else # TurnState
        return false
    end
end

# Base motile states
Run(duration::Real, speed) = RunState(; duration, speed)
Tumble(duration::Real, polar=Uniform(-π,+π),  azimuthal=Arccos()) =
    TurnState(; duration, polar, azimuthal)
Reverse(duration::Real, polar=(π,), azimuthal=Arccos()) =
    TurnState(; duration, polar, azimuthal)
Flick(duration::Real, polar=(+π/2,-π/2), azimuthal=Arccos()) =
    TurnState(; duration, polar, azimuthal)
Stop(duration::Real) =
    TurnState(; duration, polar=[zero(duration)], azimuthal=[zero(duration)])

TransitionWeights{N,T} = ProbabilityWeights{T,T,SVector{N,T}}

mutable struct Motility{N}
    current_state::Int
    motile_states::SVector{N,MotileState}
    transition_probabilities::SVector{N, TransitionWeights{N,Float64}}
end

"""
    Motility(motile_states::NTuple{N,MotileState}, rates...)
Construct a `Motility` from a sequence of `MotileState`s and their transition `rates`.

`motile_states` must be a tuple of `MotileState` instances.
`rates` must be specified in the form `(i => j, p)` where `i` is the index
of the starting state (as ordered in the `motile_states` tuple), `j`
is the index of the state towards which the transition occurs, and `p`
is the probability that, after the duration of `i` is over,
the transition occurs towards the state `j`.
It is not necessary to indicate impossible transitions; any pair
`i => j` which is not explicitly specified is assumed to have a
transition probability `p=0`.

For example, if we want to define a 3-state motility, composed by a
slow but long run, a tumble, and a fast but short run, where the two runs
can only transition towards the tumble, but the tumble can transition
with probability 30% towards the slow run and 70% towards the fast run,
we will call:

    Motility((Run(10.0, 2.0), Tumble(), Run(60.0, 0.5)),
        (1 => 2, 1.0),
        (2 => 1, 0.3), (2 => 3, 0.7),
        (3 => 2, 1.0)
    )

Some default constructors for common motility patterns are provided:
`RunTumble`, `RunReverse`, `RunReverseFlick`, `RunStop`
"""
function Motility(motile_states::NTuple{N,MotileState}, rates...) where N
    transition_matrix = zeros(N,N)
    for rate in rates
        i, j = rate[1]
        p = rate[2]
        transition_matrix[i,j] = p
    end
    transition_probabilites = SVector{N}(map(
        row -> ProbabilityWeights(SVector{N}(row)),
        eachrow(transition_matrix)
    ))
    Motility{N}(1, motile_states, transition_probabilites)
end

# Base motility patterns
function RunTumble(
    run_duration, run_speed, tumble_duration,
    polar=Uniform(-π, +π), azimuthal=Arccos()
)
    Motility(
        (Run(run_duration, run_speed), Tumble(tumble_duration, polar, azimuthal)),
        (1 => 2, 1.0), (2 => 1, 1.0)
    )
end

function RunReverse(
    run_duration_forward, run_speed_forward,
    run_duration_backward, run_speed_backward,
    reverse_duration=0.0, polar=(π,), azimuthal=Arccos()
)
    Motility(
        (
            Run(run_duration_forward, run_speed_forward),
            Reverse(reverse_duration, polar, azimuthal),
            Run(run_duration_backward, run_speed_backward),
            Reverse(reverse_duration, polar, azimuthal),
        ),
        (1 => 2, 1.0), (2 => 3, 1.0), (3 => 4, 1.0), (4 => 1, 1.0)
    )
end

function RunReverseFlick(
    run_duration_forward, run_speed_forward,
    run_duration_backward, run_speed_backward,
    reverse_duration=0.0, polar_reverse=(π,), azimuthal_reverse=Arccos(),
    flick_duration=0.0, polar_flick=(-π/2,+π/2), azimuthal_flick=Arccos(),
)
    Motility(
        (
            Run(run_duration_forward, run_speed_forward),
            Reverse(reverse_duration, polar_reverse, azimuthal_reverse),
            Run(run_duration_backward, run_speed_backward),
            Flick(flick_duration, polar_flick, azimuthal_flick)
        ),
        (1 => 2, 1.0), (2 => 3, 1.0), (3 => 4, 1.0), (4 => 1, 1.0)
    )
end

function RunStop(run_duration, run_speed, stop_duration)
    Motility(
        (Run(run_duration, run_speed), Stop(stop_duration)),
        (1 => 2, 1.0), (2 => 1, 1.0)
    )
end

# API
function update_motilestate!(microbe::AbstractMicrobe, model::AgentBasedModel)
    update_motilestate!(motilepattern(microbe), model)
    microbe.speed = rand(abmrng(model), speed(motilepattern(microbe)))
end
function update_motilestate!(motility::Motility, model::AgentBasedModel)
    i = state(motility)
    w = transition_weights(motility, i)
    j = sample(abmrng(model), eachindex(w), w)
    update_motilestate!(motility, j)
end
update_motilestate!(motility::Motility, j::Int) = (motility.current_state = j)

function motilestate(microbe::AbstractMicrobe)
    m = motilepattern(microbe)
    states(m)[state(m)]
end
state(m::Motility) = m.current_state
states(m::Motility) = m.motile_states
transition_weights(m::Motility) = m.transition_probabilities
transition_weights(m::Motility, i::Integer) = transition_weights(m)[i]
duration(m::Motility) = duration(states(m)[state(m)])
speed(m::Motility) = speed(states(m)[state(m)])
polar(m::Motility) = polar(states(m)[state(m)])
azimuthal(m::Motility) = azimuthal(states(m)[state(m)])
duration(s::MotileState) = s.duration
speed(s::MotileState) = s.speed
polar(s::MotileState) = s.polar
azimuthal(s::MotileState) = s.azimuthal
