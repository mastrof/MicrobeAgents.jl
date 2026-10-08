export Behavior, behaviors, findbehavior
export initialize!, affect!, bias, speed_factor, transition_weights!

"""
    initialize!(behavior, microbe, model)
Hook called once when `microbe` is added to `model`. Default: no-op.
"""
initialize!(b, microbe, model) = nothing

"""
    affect!(behavior, microbe, model)
Hook called at every step right after translation, to update internal state
or apply arbitrary per-step logic. Default: no-op.
A plain function `f` used as a behavior is called as `f(microbe, model)`.
"""
affect!(b, microbe, model) = nothing
affect!(f::Function, microbe, model) = (f(microbe, model); nothing)

"""
    bias(behavior, microbe, model)
Hook returning a multiplicative factor on the rate of switching out of
biased motile states (see `biased`). Default: `1.0`.

    bias(microbe, model)
Total bias of `microbe`, i.e. the product of `bias` over its behaviors.
"""
bias(b, microbe, model) = 1.0
bias(microbe::AbstractMicrobe, model) =
    _prod(b -> bias(b, microbe, model), values(behaviors(microbe)))

"""
    speed_factor(behavior, microbe)
Hook returning a multiplicative factor on the microbe speed. Default: `1.0`.
Must only depend on state stored in `behavior` or `microbe`.
"""
speed_factor(b, microbe) = 1.0

"""
    transition_weights!(w, behavior, microbe, model)
Hook to modify `w`, a scratch copy of the transition weights out of the current
motile state, right before the next state is sampled. Default: no-op.
"""
transition_weights!(w, b, microbe, model) = nothing

"""
    Behavior(; affect!, bias, speed_factor, transition_weights!)
Behavior built from functions, for stateless behaviors without defining a type.
Each keyword is optional; signatures are those of the hooks without the
behavior argument, e.g. `bias = (microbe, model) -> 2.0`.
"""
struct Behavior{A,B,S,T}
    affect::A
    bias::B
    speed_factor::S
    transition_weights::T
end
Behavior(; affect! = nothing, bias = nothing, speed_factor = nothing,
    transition_weights! = nothing) =
    Behavior(affect!, bias, speed_factor, transition_weights!)

_call(::Nothing, default, args...) = default
_call(f, default, args...) = f(args...)
affect!(b::Behavior, microbe, model) = (_call(b.affect, nothing, microbe, model); nothing)
bias(b::Behavior, microbe, model) = _call(b.bias, 1.0, microbe, model)
speed_factor(b::Behavior, microbe) = _call(b.speed_factor, 1.0, microbe)
transition_weights!(w, b::Behavior, microbe, model) =
    (_call(b.transition_weights, nothing, w, microbe, model); nothing)

"""
    behaviors(microbe)
Return the behaviors (`Tuple` or `NamedTuple`) of `microbe`.
"""
behaviors(m::AbstractMicrobe) = ()

"""
    findbehavior(microbe, T)
Return the first behavior of `microbe` of type `T`, or `nothing`.
"""
findbehavior(m::AbstractMicrobe, ::Type{T}) where {T} = _findfirst(T, values(behaviors(m)))
_findfirst(::Type, ::Tuple{}) = nothing
_findfirst(::Type{T}, t::Tuple) where {T} =
    first(t) isa T ? first(t) : _findfirst(T, Base.tail(t))

# compile-time unrolled loops over heterogeneous tuples
@inline _foreach(f, ::Tuple{}) = nothing
@inline _foreach(f, t::Tuple) = (f(first(t)); _foreach(f, Base.tail(t)))
@inline _prod(f, ::Tuple{}) = 1.0
@inline _prod(f, t::Tuple) = f(first(t)) * _prod(f, Base.tail(t))

# each agent gets its own copy of stateful behaviors;
# functions and `Behavior`s are stateless and shared
copybehavior(b) = deepcopy(b)
copybehavior(f::Function) = f
copybehavior(b::Behavior) = b
copybehaviors(bs::Union{Tuple,NamedTuple}) = map(copybehavior, bs)
copybehaviors(x) = throw(ArgumentError(
    "`behaviors` must be a Tuple or NamedTuple, got $(typeof(x)); " *
    "wrap single behaviors as `behaviors = (b,)`"
))

initialize_behaviors!(microbe, model) =
    _foreach(b -> initialize!(b, microbe, model), values(behaviors(microbe)))
