# Composable microbe behaviors — design

Date: 2026-10-08
Status: approved in brainstorming, pending spec review

## Goal

Make motile and chemotactic behaviors freely combinable and modifiable without
defining new agent types. Concretely:

- combine sensing models with motility modulations (chemokinesis, speed-dependent
  turn rates, state-dependent transition weights);
- add non-chemical behaviors (walls/surfaces, other cues, interactions);
- ultimately, plug in arbitrary user-defined per-step logic.

Constraints from the user:

- breaking changes are acceptable, but model creation (`StandardABM(...)`) and agent
  creation (`add_agent!(...)`) should change as little as possible;
- the user-facing API must not require declaring or managing types. Users call
  constructors and write functions; struct definitions are only needed for
  behaviors with per-agent state.

## Problem with the current design

Behavior and parameters are bound to the agent type. Each chemotaxis model is an
`@agent struct` whose fields hold its parameters, and its behavior is selected by
dispatching `chemotaxis!`/`bias` on that type. Adding one parameter or one rule
means a new agent type. `SonMenolascina` shows the cost: it is four behaviors
(BrownBerg sensing, chemokinesis, speed-dependent bias, speed-dependent flick
probability) hard-wired into one type via overloads of `speed`,
`switching_probability` and `update_motilestate!`, none reusable on their own.
`model.affect!` is model-global, so per-agent behavioral differences require
different types (and `Union` model types).

## Design overview

A microbe holds a tuple (or NamedTuple) of **behaviors**. A behavior is any value;
what it does is defined by which **hooks** are implemented for its type. The
stepping pipeline calls each hook on every behavior, in tuple order (loop unrolled
at compile time). Names in a NamedTuple are arbitrary labels used only for access.

```julia
model = StandardABM(Microbe{2}, space, dt; properties)   # unchanged
add_agent!(model;
    motility = RunReverseFlick(...),
    rotational_diffusivity = 0.035, radius = 0.5,
    behaviors = (
        chemotaxis = BrownBerg(gain = 600),
        kinesis    = Chemokinesis(threshold = 0.05, factor = 1.3),
        walls      = my_wall_rule!,                           # plain function
        light      = Behavior(bias = (m, model) -> 1 + 2*light(position(m))),
    ),
)
```

## Hook protocol

| Hook | Called | Default | Combination |
|---|---|---|---|
| `initialize!(b, microbe, model)` | once in `add_agent!`, after placement | no-op | tuple order |
| `affect!(b, microbe, model)` | every step, right after move | no-op | tuple order |
| `bias(b, microbe, model)` | every step, only if `biased(motilestate(microbe))` | `1.0` | product |
| `speed_factor(b, microbe)` | whenever `speed(microbe)` is read | `1.0` | product |
| `transition_weights!(w, b, microbe, model)` | only when a switch occurs, before sampling | no-op | tuple order |

Rules:

- A plain `Function` used as a behavior is an `affect!`: called as `f(microbe, model)`.
- `Behavior(; affect!, bias, speed_factor, transition_weights!)` wraps functions
  (each optional, signatures as above minus the leading `b`) so stateless
  behaviors with any hook need no type declaration. It holds no per-agent state;
  behaviors needing memory between steps are user structs (typically
  `@kwdef mutable struct`) with methods for the needed hooks.
- `speed_factor` has no `model` argument because `speed`/`velocity` are called in
  contexts without the model (including inside sensing). Its value must come from
  the behavior's stored state. `microbe.speed` remains the sampled base speed;
  `speed(microbe) = microbe.speed * prod(speed_factor(b, microbe) for b in behaviors)`.
- `transition_weights!` receives `w`, a scratch copy of the current state's
  transition weights; the `Motility` itself is never mutated by behaviors.
- Aggregate accessor: `bias(microbe, model)` returns the product over behaviors
  (used by `switching_probability` and for data collection).

## Stepping pipeline

```
microbe_step!(microbe, model):
    move_step!(microbe, model)          # unchanged: translation + rotational diffusion
    affect_step!(microbe, model)        # invalidate field cache; foreach affect!
    reorient_step!(microbe, model):
        if can_turn(microbe): turn!     # unchanged
        p = switching_probability       # = (biased ? bias(microbe, model) : 1) * dt / τ
        if rand < p:
            w = copy(base weights of current state)
            foreach transition_weights!(w, b, microbe, model)
            sample next state from w
            if new state is a zero-duration TurnState:
                turn!; repeat the weighted sampling above (hooks applied again,
                consistent with today's overridden update_motilestate!)
            update_speed!
```

`microbe_pathfinder_step!` is the same with `move_step_pathfinder!`.
No pre-move hook for now; `affect!` runs immediately after translation, which
covers boundary/wall rules. Replacing translation remains an `agent_step!` choice.

## Memoized field quantities

Environmental quantities needed by several behaviors are evaluated at most once
per microbe step, lazily:

- the model holds one `FieldCache{D}` (internal property) with a slot per built-in
  quantity: `concentration::Union{Nothing,Float64}`,
  `gradient::Union{Nothing,SVector{D,Float64}}`, `time_derivative`,
  `chemoattractant_diffusivity`; plus the id of the microbe it refers to and a
  validity flag;
- `affect_step!` resets the cache for the current microbe; it is marked invalid at
  the end of `microbe_step!`;
- `concentration(microbe, model)` & co. (existing microbe-level wrappers in
  `fields.jl`) read the cache when valid for this microbe, otherwise compute and
  store. Calls outside a step (e.g. from `adata`) always compute fresh values;
- the `AbstractChemoattractant` interface (`f(microbe, model)` functions) is unchanged;
- the cache is reset in `affect_step!` and invalidated in `move_step!` and at the
  end of `reorient_step!`, so reordering subroutines never yields stale values;
- custom expensive quantities (obstacle distance, neighbour searches, extra
  fields) are handled by explicit *sensor* behaviors placed early in the tuple,
  storing the value in their own field; other behaviors read it via the
  NamedTuple name or `findbehavior(microbe, T)` (first behavior of type `T`, or
  `nothing`). No generic `Dict`-based cache (type-unstable, allocating).

Assumes sequential stepping per model (as Agents.jl `StandardABM` does).

## Agent type and creation

- `Microbe{D,N,B} <: AbstractMicrobe{D,N}` fields: `id, pos, vel, speed`,
  `motility::Motility{N}`, `rotational_diffusivity = 0.0`, `radius = 0.0`,
  `behaviors::B = ()`. The generic `state` field is removed.
- `AbstractMicrobe{D,N}` stays; the pipeline uses the accessor `behaviors(m)`
  (default `()`), so custom agent types keep working.
- `StandardABM(Microbe{D}, space, dt; ...)` unchanged; the `:affect!` default
  property is removed, the field cache is added. Mixed populations need no
  `Union`: all microbes are `Microbe{D,…}`.
- `add_agent!(...; motility, behaviors = (), kw...)`: `B` inferred like `N`;
  stateful behaviors are `deepcopy`'d per agent so state is never shared
  (plain functions and `Behavior` wrappers are shared, not copied); a non-Tuple
  `behaviors` value raises an `ArgumentError`; `initialize!` is called for each
  behavior right before the agent is placed in the model, so a throwing
  `initialize!` leaves the model untouched.

## Library behaviors (migration of existing models)

The chemotaxis agent types become behavior structs with the same names, holding
only sensing parameters and internal state:

| Behavior | Fields (defaults as today) | Hooks |
|---|---|---|
| `BrownBerg` | `gain=660, receptor_binding_constant=100, memory=1, state=0` | `affect!`, `bias` |
| `Brumley` | `memory=1.3, gain_receptor=50, gain=50, chemotactic_precision=6, state=0` | `initialize!` (radius check), `affect!`, `bias` |
| `Celani` | `gain=50, memory=1, chemotactic_precision=0, markovian_variables, state=0` | `initialize!` (steady state + radius check), `affect!`, `bias` |
| `Xie` | current parameters, `state, state_m, state_z` | `initialize!` (radius check), `affect!`, `bias` (forward/backward by motile state index, as today) |
| `Chemokinesis` | `threshold=0.05, factor=1.3, on=false` | `affect!`, `speed_factor` |
| `SpeedDependentTurnRate` | `eta, ζ, θ, vT` (SonMenolascina values) | `bias` |
| `SpeedDependentFlick` | logistic parameters currently hard-coded | `transition_weights!` (only for 4-state motilities, when leaving the backward run, state 3; no-op otherwise, as today) |

- `SonMenolascina(; kw...)` becomes a function returning the NamedTuple
  `(chemokinesis = Chemokinesis(...), chemotaxis = BrownBerg(...),
  turnrate = SpeedDependentTurnRate(...), flick = SpeedDependentFlick())`
  with SonMenolascina's parameter values; usable as `behaviors = SonMenolascina()`
  or splatted into a larger tuple. `Chemokinesis` comes first because the original
  model updated chemokinesis before BrownBerg read the (speed-dependent) velocity.
- Xie's unused `turn_rate_forward`/`turn_rate_backward` fields are dropped.
- Physical properties (`rotational_diffusivity`, `radius`) and motility are passed
  to `add_agent!`; per-model defaults (e.g. BrownBerg's 0.035 rad²/s) are dropped.
- `initialize!` of noisy models throws an informative error if
  `chemotactic_precision > 0` and `radius(microbe) == 0` (noise would be `Inf`).
- Celani's custom `add_agent!` methods are removed (replaced by `initialize!`).

## Removed / changed public API

Removed: `chemotaxis!`, `bias(microbe)` single-argument methods on agent types,
`state(microbe)` for microbes, the `:affect!` model property, the chemotaxis
agent types (names reused for behaviors), `SonMenolascina` agent type.

Added (exported, with docstrings): `Behavior`, `behaviors`, `findbehavior`,
`initialize!`, `affect!`, `bias` (3-arg hook and 2-arg aggregate), `speed_factor`,
`transition_weights!`, `Chemokinesis`, `SpeedDependentTurnRate`,
`SpeedDependentFlick`.

Unchanged: `microbe_step!`, `microbe_pathfinder_step!`, `move_step!`,
`move_step_pathfinder!`, `affect_step!`, `reorient_step!`,
`switching_probability`, `can_turn`, `Motility` and constructors, fields API,
Analysis submodule.

Data collection: internal state via `adata = [m -> m.behaviors.chemotaxis.state]`
(or positional index); aggregate bias via `adata = [m -> bias(m, model)]`.

## Testing

- Hook unit tests: neutral defaults; product combination of `bias` and
  `speed_factor`; `transition_weights!` operates on a scratch copy (motility
  weights unchanged after a switch); per-agent behavior copies are independent;
  plain functions and `Behavior` wrappers dispatch correctly; `findbehavior`.
- Cache tests: a counting chemoattractant shows one evaluation per quantity per
  microbe step when several behaviors query it; fresh evaluation outside steps.
- Regression: before refactoring, record reference trajectories/state/bias time
  series on `main` for each chemotaxis model (fixed seed, small model) and store
  them as test fixtures; migrated behaviors must reproduce them exactly (RNG draw
  order is preserved by construction).
- Existing tests migrated to the new API; docs build runs all examples.

## Open details (to settle in the implementation plan)

- `Base.show` for `Microbe` should list behaviors by name/type.
