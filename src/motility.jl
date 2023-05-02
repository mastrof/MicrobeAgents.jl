export
    AbstractMotility, AbstractMotilityOneStep, AbstractMotilityTwoStep,
    MotilityOneStep, MotilityTwoStep,
    RunTumble, RunReverse, RunReverseFlick,
    MotileState, TwoState, Forward, Backward, switch!,
    rand_speed, rand_polar, rand_azimuthal

"""
    AbstractMotility
General abstract interface for motility patterns.
"""
abstract type AbstractMotility end



# Copied and adapted from the @agent macro of Agents.jl
macro motility(new_name, base_type, extra_fields)
    # This macro was generated with the guidance of @rdeits on Discourse:
    # https://discourse.julialang.org/t/
    # metaprogramming-obtain-actual-type-from-symbol-for-field-inheritance/84912

    # We start with a quote. All macros return a quote to be evaluated
    quote
        let
            # Here we collect the field names and types from the base type
            # Because the base type already exists, we escape the symbols to obtain it
            base_fieldnames = fieldnames($(esc(base_type)))
            base_fieldtypes = [t for t in getproperty($(esc(base_type)), :types)]
            base_fields = [:($f::$T) for (f, T) in zip(base_fieldnames, base_fieldtypes)]
            # Then, we prime the additional name and fields into QuoteNodes
            # We have to do this to be able to interpolate them into an inner quote.
            name = $(QuoteNode(new_name))
            additional_fields = $(QuoteNode(extra_fields.args))
            # Now we start an inner quote. This is because our macro needs to call `eval`
            # However, this should never happen inside the main body of a macro
            # There are several reasons for that, see the cited discussion at the top
            expr = quote
                mutable struct $name <: AbstractMotility
                    $(base_fields...)
                    $(additional_fields...)
                end
            end
            # @show expr # uncomment this to see that the final expression looks as desired
            # It is important to evaluate the macro in the module that it was called at
            Core.eval($(__module__), expr)
        end
        # allow attaching docstrings to the new struct, issue #715
        Core.@__doc__($(esc(Docs.namify(new_name))))
        nothing
    end
end
# There should be away that only the 4-argument version is used
# and the 3-argument version just passes `AbstractAgent` to the 4-argument.
macro motility(new_name, base_type, super_type, extra_fields)
    # This macro was generated with the guidance of @rdeits on Discourse:
    # https://discourse.julialang.org/t/
    # metaprogramming-obtain-actual-type-from-symbol-for-field-inheritance/84912

    # We start with a quote. All macros return a quote to be evaluated
    quote
        let
            # Here we collect the field names and types from the base type
            # Because the base type already exists, we escape the symbols to obtain it
            base_fieldnames = fieldnames($(esc(base_type)))
            base_fieldtypes = [t for t in getproperty($(esc(base_type)), :types)]
            base_fields = [:($f::$T) for (f, T) in zip(base_fieldnames, base_fieldtypes)]
            # Then, we prime the additional name and fields into QuoteNodes
            # We have to do this to be able to interpolate them into an inner quote.
            name = $(QuoteNode(new_name))
            additional_fields = $(QuoteNode(extra_fields.args))
            # Now we start an inner quote. This is because our macro needs to call `eval`
            # However, this should never happen inside the main body of a macro
            # There are several reasons for that, see the cited discussion at the top
            expr = quote
                # Also notice that we escape supertype and interpolate it twice
                # because this is expected to already be defined in the calling module
                mutable struct $name <: $$(esc(super_type))
                    $(base_fields...)
                    $(additional_fields...)
                end
            end
            # @show expr # uncomment this to see that the final expression looks as desired
            # It is important to evaluate the macro in the module that it was called at
            Core.eval($(__module__), expr)
        end
        # allow attaching docstrings to the new struct, issue #715
        Core.@__doc__($(esc(Docs.namify(new_name))))
        nothing
    end
end




"""
    AbstractMotilityOneStep
One-step motility patterns (`RunTumble`).
Subtypes have the following fields:
- `speed`: distribution of microbe speed, new values extracted after each turn
- `polar`: distribution of polar angles
- `azimuthal`: distribution azimuthal angles
For 2-dimensional microbe types, only `polar` defines reorientations and `azimuthal` is ignored.

New one-step motility patterns can be created as
```
MicrobeAgents.@motility NewMotilityType MotilityOneStep AbstractMotilityOneStep begin
    # some extra fields if needed
end
```
The necessary fields and subtyping will be added automatically, only new extra fields
need to be specified in the definition.
If default values have to be specified, a constructor needs to be defined explicitly.
"""
abstract type AbstractMotilityOneStep <: AbstractMotility end
# only used to define the base fields
struct MotilityOneStep
    speed
    polar
    azimuthal
end

@motility RunTumble MotilityOneStep AbstractMotilityOneStep begin; end
"""
    RunTumble(; speed=(30.0,), polar=Uniform(-π,π), azimuthal=Arccos())
Run-tumble motility.
The kwargs `speed`, `polar` and `azimuthal` must be sampleable objects
(ranges, arrays, tuples, distributions...), _not scalars_.

With default values, the reorientations are uniform on the sphere.
"""
RunTumble(; speed=(30.0,), polar=Uniform(-π,π), azimuthal=Arccos()) = RunTumble(speed, polar, azimuthal)


"""
    AbstractMotilityTwoStep
Two-step motility patterns (`RunReverse` and `RunReverseFlick`), with different
properties between forward and backward state of motion.
Subtypes have the following fields:
- `speed_forward`: distribution of microbe speed, new values extracted after each turn
- `polar_forward`: distribution of in-plane reorientations for motile state fw
- `azimuthal_forward`: distribution of out-of-plane reorientations for motile state fw
- `speed_backward`: distribution of microbe speed, new values extracted after each turn
- `polar_backward`: distribution of in-plane reorientations for motile state bw
- `azimuthal_backward`: distribution of out-of-plane reorientations for motile state bw
- `motile_state`: defines current motile state (e.g. `Forward` or `Backward` for a `TwoState`)
For 2-dimensional microbe types, only polar distributions define reorientations
while azimuthal ones are ignored.

New two-step motility patterns can be created as
```
MicrobeAgents.@motility NewMotilityType MotilityTwoStep AbstractMotilityTwoStep begin
    # some extra fields if needed
end
```
The necessary fields and subtyping will be added automatically, only new extra fields
need to be specified in the definition.
If default values have to be specified, a constructor needs to be defined explicitly.
"""
abstract type AbstractMotilityTwoStep <: AbstractMotility end
# only used to define the base fields
struct MotilityTwoStep
    speed_forward
    polar_forward
    azimuthal_forward
    speed_backward
    polar_backward
    azimuthal_backward
    motile_state
end

@motility RunReverse MotilityTwoStep AbstractMotilityTwoStep begin; end
"""
    RunReverse(;
        speed_forward = (30.0,),
        polar_forward = (π,),
        azimuthal_forward = Arccos(),
        speed_backward = speed_forward,
        polar_backward = polar_forward,
        azimuthal_backward = azimuthal_forward
    )
Run-reverse motility, with possibility to have different properties between
the forward (run) and backward (reverse) stages.
All the fields must be sampleable objects (ranges, arrays, tuples, distributions...),
_not scalars_.

With default values, reorientations are always perfect reversals
and the speed is identical between forward and backward runs.
"""
RunReverse(;
    speed_forward = (30.0,),
    polar_forward = (π,),
    azimuthal_forward = Arccos(),
    speed_backward = speed_forward,
    polar_backward = polar_forward,
    azimuthal_backward = azimuthal_forward,
    motile_state = MotileState()
) = RunReverse(
    speed_forward, polar_forward, azimuthal_forward,
    speed_backward, polar_backward, azimuthal_backward,
    motile_state)


@motility RunReverseFlick MotilityTwoStep AbstractMotilityTwoStep begin; end
"""
    RunReverseFlick(;
        speed_forward = (30.0,),
        polar_forward = (π,),
        azimuthal_forward = Arccos(),
        speed_backward = speed_forward,
        polar_backward = (-π/2, π/2),
        azimuthal_backward = azimuthal_forward
    )
Run-reverse-flick motility, with possibility to have different properties between
the forward (run) and backward (reverse) stages.
All the fields must be sampleable objects (ranges, arrays, tuples, distributions...),
_not scalars_.

With default values, reorientations after forward runs are perfect reversals,
while reorientations after backward runs are uniformly distributed on the circle
normal to the run direction; speed is identical between forward and backward runs.
"""
RunReverseFlick(;
    speed_forward = (30.0,),
    polar_forward = (π,),
    azimuthal_forward = Arccos(),
    speed_backward = speed_forward,
    polar_backward = (-π/2, π/2),
    azimuthal_backward = azimuthal_forward,
    motile_state = MotileState()
) = RunReverseFlick(
    speed_forward, polar_forward, azimuthal_forward,
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
function Base.getproperty(obj::AbstractMotilityTwoStep, sym::Symbol)
    if sym === :state
        return obj.motile_state.state
    else
        return getfield(obj, sym)
    end
end
function Base.setproperty!(obj::AbstractMotilityTwoStep, sym::Symbol, x)
    if sym === :state
        return setfield!(obj.motile_state, :state, x)
    else
        return setfield!(obj, sym, x)
    end
end
# define rules for switching motile state
switch!(::AbstractMotilityOneStep) = nothing
""" 
    switch!(m::AbstractMotilityTwoStep)
Switch the state of a two-step motility pattern (`m.state`)
from `Forward` to `Backward` and viceversa.
"""
switch!(m::AbstractMotilityTwoStep) = (m.state = ~m.state; nothing)
Base.:~(x::TwoState) = TwoState(~Bool(x))

# convenient wrapper to sample new motile state
Base.rand(m::AbstractMotility) = (rand_speed(m), rand_polar(m), rand_azimuthal(m))
Base.rand(rng::AbstractRNG, m::AbstractMotility) = (rand_speed(rng,m), rand_polar(rng,m), rand_azimuthal(rng,m))
# expand rand_vel from utils.jl
"""
    rand_speed([rng,] m::AbstractMotilityOneStep)
Extract value from the speed distribution of the motility pattern `m.speed`.
"""
rand_speed(m::AbstractMotilityOneStep) = rand(m.speed)
rand_speed(rng::AbstractRNG, m::AbstractMotilityOneStep) = rand(rng, m.speed)
"""
    rand_speed([rng,] m::AbstractMotilityTwoStep)
Extract value from the speed distribution of the motility pattern.
If `motilestate(m) == ForwardState()` extract from `speed_forward`, otherwise
from `speed_backward`.
"""
function rand_speed(m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(m.speed_forward)
    else
        return rand(m.speed_backward)
    end
end
function rand_speed(rng::AbstractRNG, m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(rng, m.speed_forward)
    else
        return rand(rng, m.speed_backward)
    end
end
"""
    rand_vel([rng,] N::Int, m::AbstractMotility)
Generate a random N-tuple, with norm defined by the speed distribution of
the motile pattern `m`.
"""
rand_vel(D::Int, m::AbstractMotility) = rand_vel(D) .* rand_speed(m)
rand_vel(rng::AbstractRNG, D::Int, m::AbstractMotility) = rand_vel(rng, D) .* rand_speed(m)

# extract angles from polar and azimuthal distributions
rand_polar(m::AbstractMotilityOneStep) = rand(m.polar)
rand_polar(rng::AbstractRNG, m::AbstractMotilityOneStep) = rand(rng, m.polar)
rand_azimuthal(m::AbstractMotilityOneStep) = rand(m.azimuthal)
rand_azimuthal(rng::AbstractRNG, m::AbstractMotilityOneStep) = rand(rng, m.azimuthal)
function rand_polar(m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(m.polar_forward)
    else
        return rand(m.polar_backward)
    end
end
function rand_polar(rng::AbstractRNG, m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(rng, m.polar_forward)
    else
        return rand(rng, m.polar_backward)
    end
end
function rand_azimuthal(m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(m.azimuthal_forward)
    else
        return rand(m.azimuthal_backward)
    end
end
function rand_azimuthal(rng::AbstractRNG, m::AbstractMotilityTwoStep)
    if m.state == Forward
        return rand(rng, m.azimuthal_forward)
    else
        return rand(rng, m.azimuthal_backward)
    end
end