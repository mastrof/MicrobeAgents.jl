export MotileState, RunState, TurnState, Run, Tumble, Reverse, Flick, Stop
export Motility, RunTumble, RunReverse, RunReverseFlick, RunStop
export update_motilestate!
export motilepattern, motilestate, state, states, transition_weights
export duration, speed, polar, azimuthal
export Arccos # from Agents

@sum_structs MotileState begin
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

@pattern biased(::RunState) = true
@pattern biased(::TurnState) = false

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
        row -> ProbabilityWeights(MVector{N}(row)),
        eachrow(transition_matrix)
    ))
    Motility{N}(1, motile_states, transition_probabilites)
end

# Base motility patterns
function RunTumble(
    run_duration, run_speed,
    tumble_duration=0, polar=Uniform(-π, +π), azimuthal=Arccos()
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
    reverse_duration=0, flick_duration=0,
    polar_reverse=(π,), azimuthal_reverse=Arccos(),
    polar_flick=(-π/2,+π/2), azimuthal_flick=Arccos(),
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
    polar(motility::Motility)
Return the distribution of polar reorientation angles of the
current motile state.
"""
polar(m::Motility) = polar(states(m)[state(m)])
"""
    azimuthal(motility::Motility)
Return the distribution of azimuthal reorientation angles of the
current motile state.
"""
azimuthal(m::Motility) = azimuthal(states(m)[state(m)])
duration(s::MotileState) = s.duration
speed(s::MotileState) = s.speed
polar(s::MotileState) = s.polar
azimuthal(s::MotileState) = s.azimuthal
