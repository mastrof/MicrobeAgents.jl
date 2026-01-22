export MotileState, RunState, TurnState
export biased
export Run, Tumble, Reverse, Flick, Stop
export Motility, RunTumble, RunReverse, RunReverseFlick, RunStop
export update_motilestate!, update_speed!
export motilepattern, motilestate, state, states, transition_weights
export duration, speed, angle
export Arccos # from Agents

struct RunState
    duration
    speed
    angle
    azimuthal
end
struct TurnState
    duration
    speed
    angle
    azimuthal
end
@sumtype MotileState(RunState, TurnState)
function RunState(; duration, speed, angle=nothing, azimuthal=nothing)
    MotileState(RunState(duration, speed, angle, azimuthal))
end
function TurnState(; duration, angle, speed=[zero(duration)], azimuthal=Uniform(-π,π))
    MotileState(TurnState(duration, speed, angle, azimuthal))
end

"""
    biased(s::MotileState)
Check whether a motile state can have its duration biased
by an internal state.

MicrobeAgents implements two motile state types:
`RunState` and `TurnState`.
`RunState` is biased (i.e., this function returns `true`),
`TurnState` is not.
"""
biased(s::MotileState) = biased(variant(s))
biased(::RunState) = true
biased(::TurnState) = false

# Base motile states
Run(speed, duration) = RunState(; duration, speed)
Tumble(angle, duration=0.0) = TurnState(; duration, angle)
Reverse(angle=(π,), duration=0.0) = TurnState(; duration, angle)
Flick(angle=(-π/2, π/2), duration=0.0) = TurnState(; duration, angle)
Stop(duration) =
    TurnState(; duration, angle=[zero(duration)], azimuthal=[zero(duration)])

TransitionWeights{N,T} = ProbabilityWeights{T,T,MVector{N,T}}

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
slow but long run, an instantaneous tumble with some `angle_pdf`,
and a fast but short run, where the two runs
can only transition towards the tumble, but the tumble can transition
with probability 30% towards the slow run and 70% towards the fast run,
we will call:

    Motility((Run(10.0, 2.0), Tumble(0.0, angle_pdf), Run(60.0, 0.5)),
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
        row -> ProbabilityWeights(MVector{N}(row)),
        eachrow(transition_matrix)
    ))
    Motility{N}(1, motile_states, transition_probabilites)
end

# Base motility patterns
"""
    RunTumble(run_speed, run_duration, angle; tumble_duration=0.0)
"""
RunTumble(
    run_speed, run_duration, angle; kwargs...
) = RunTumble(;
    run_speed, run_duration, angle, kwargs...
)
function RunTumble(;
    run_speed, run_duration,
    angle, tumble_duration=0.0,
)
    Motility(
        (Run(run_speed, run_duration), Tumble(angle, tumble_duration)),
        (1 => 2, 1.0), (2 => 1, 1.0)
    )
end

"""
    RunReverse(
        run_speed_forward, run_duration_forward,
        run_speed_backward, run_duration_backward;
        angle=(π,), reverse_duration=0.0,
"""
RunReverse(
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward,
    kwargs...
) = RunReverse(;
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward,
    kwargs...
)
function RunReverse(;
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward,
    angle=(π,), reverse_duration=0.0,
)
    Motility(
        (
            Run(run_speed_forward, run_duration_forward),
            Reverse(angle, reverse_duration),
            Run(run_speed_backward, run_duration_backward),
            Reverse(angle, reverse_duration),
        ),
        (1 => 2, 1.0), (2 => 3, 1.0), (3 => 4, 1.0), (4 => 1, 1.0)
    )
end

"""
    RunReverseFlick(
        run_speed_forward, run_duration_forward,
        run_speed_backward, run_duration_backward;
        reverse_angle=(π,) flick_angle=(-π/2, +π/2),
        reverse_duration=0.0, flick_duration=0.0,
    )
"""
RunReverseFlick(
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward;
    kwargs...
) = RunReverseFlick(;
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward,
    kwargs...
)
function RunReverseFlick(;
    run_speed_forward, run_duration_forward,
    run_speed_backward, run_duration_backward,
    reverse_angle=(π,), flick_angle=(-π/2,+π/2),
    reverse_duration=0, flick_duration=0,
)
    Motility(
        (
            Run(run_speed_forward, run_duration_forward),
            Reverse(reverse_angle, reverse_duration),
            Run(run_speed_backward, run_duration_backward),
            Flick(flick_angle, flick_duration)
        ),
        (1 => 2, 1.0), (2 => 3, 1.0), (3 => 4, 1.0), (4 => 1, 1.0)
    )
end

"""
    RunStop(run_speed, run_duration, stop_duration)
"""
RunStop(
    run_speed, run_duration, stop_duration
) = RunStop(;
    run_speed, run_duration, stop_duration
)
function RunStop(; run_speed, run_duration, stop_duration)
    Motility(
        (Run(run_speed, run_duration), Stop(stop_duration)),
        (1 => 2, 1.0), (2 => 1, 1.0)
    )
end

# API
"""
    update_motilestate!(microbe, model)
Update the motile state of `microbe` by randomly sampling the next state
according to the prescribed transition weights.
"""
function update_motilestate!(microbe::AbstractMicrobe, model::AgentBasedModel)
    update_motilestate!(motilepattern(microbe), model)
end
function update_motilestate!(motility::Motility, model::AgentBasedModel)
    i = state(motility)
    w = transition_weights(motility, i)
    j = sample(abmrng(model), eachindex(w), w)
    update_motilestate!(motility, j)
end
update_motilestate!(motility::Motility, j::Int) = (motility.current_state = j)

"""
    update_speed!(microbe, model)
Update the speed of `microbe` by randomly sampling from the
speed distribution of the current motile state.
"""
function update_speed!(microbe::AbstractMicrobe, model::AgentBasedModel)
    microbe.speed = rand(abmrng(model), speed(motilestate(microbe)))
end

"""
    motilestate(microbe)
Return the current motile state of `microbe`.
"""
function motilestate(microbe::AbstractMicrobe)
    m = motilepattern(microbe)
    states(m)[state(m)]
end
"""
    state(motility::Motility)
Return the index of active motile state.
"""
state(m::Motility) = m.current_state
"""
    state(motility::Motility)
Return the list of all motile states.
"""
states(m::Motility) = m.motile_states
transition_weights(m::Motility) = m.transition_probabilities
"""
    transition_weights(motility::Motility, i::Integer)
Return the transition weights from the state with index `i`
towards the other motile states.
"""
transition_weights(m::Motility, i::Integer) = transition_weights(m)[i]
"""
    duration(motility::Motility)
Return the average unbiased duration of the current motile state.
"""
duration(m::Motility) = duration(states(m)[state(m)])
"""
    speed(motility::Motility)
Return the speed distribution of the current motile state.
"""
speed(m::Motility) = speed(states(m)[state(m)])
"""
    angle(motility::Motility)
Return the distribution of polar reorientation angles of the
current motile state.
"""
angle(m::Motility) = angle(states(m)[state(m)])
"""
    azimuthal(motility::Motility)
Return the distribution of azimuthal reorientation angles of the
current motile state.
"""
azimuthal(m::Motility) = azimuthal(states(m)[state(m)])
duration(s::MotileState) = s.duration
speed(s::MotileState) = s.speed
angle(s::MotileState) = s.angle
azimuthal(s::MotileState) = s.azimuthal
