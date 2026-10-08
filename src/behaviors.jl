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
A plain `Function` `f` used as a behavior is called as `f(microbe, model)`;
a callable struct needs its own `affect!` method.
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
Modify `w` by indexing (`w[i] = x`), which keeps `sum(w)` consistent; do not
edit `w.values` directly.
"""
transition_weights!(w, b, microbe, model) = nothing

"""
    Behavior(; affect!, bias, speed_factor, transition_weights!)
Behavior built from functions, for stateless behaviors without defining a type.
Each keyword is optional; signatures are those of the hooks without the
behavior argument, e.g. `bias = (microbe, model) -> 2.0`.

`add_agent!` copies behaviors for each agent: struct behaviors are `deepcopy`'d
(including any data they reference), while plain functions and `Behavior`
objects are shared, so mutable state captured by their closures is shared too.
To share data between agents, store it in the model properties, or opt out of
copying for a type with `MicrobeAgents.copybehavior(b::MyType) = b`.
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

"""
    MicrobeAgents.copybehavior(behavior)
Extension point: how a behavior is copied for each new agent.
Default: `deepcopy` (functions and `Behavior`s are returned as is).
Define `MicrobeAgents.copybehavior(b::MyType) = b` to share one instance.
"""
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
