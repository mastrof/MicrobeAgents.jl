# Behaviors

Everything a microbe does beyond its motility pattern is expressed by
**behaviors**: values passed to `add_agent!` through the `behaviors` keyword,
as a `Tuple` or a `NamedTuple`. Names are arbitrary labels, useful to access
the behaviors later; the order of the behaviors is the order of execution.

```julia
model = StandardABM(Microbe{2}, space, dt; properties)
add_agent!(model;
    motility = RunReverseFlick([30.0], 0.5, [30.0], 0.5),
    radius = 0.5,
    behaviors = (
        chemotaxis = BrownBerg(gain = 600),
        kinesis = Chemokinesis(threshold = 0.05, factor = 1.3),
    ),
)
model[1].behaviors.chemotaxis.state # internal state of the chemotaxis behavior
```

Each agent receives its own copy of the behaviors. Struct behaviors are
`deepcopy`'d, including any data they reference (large arrays, obstacle lists,
shared fields): to share such data, reference it from the model properties, or
opt out of copying for a type with `MicrobeAgents.copybehavior(b::MyType) = b`.
Plain functions and `Behavior` objects are shared across agents, so any mutable
state captured by their closures is shared too.

## Hooks

During each step the microbe calls these hooks on each behavior.
A behavior only defines the hooks it needs; the others are neutral.

```@docs
initialize!
affect!
bias
speed_factor
transition_weights!
```

The step is: translation and rotational diffusion (`move_step!`), then
`affect!` for each behavior (`affect_step!`), then active reorientation and
motile state switching (`reorient_step!`), where the switching rate is
multiplied by the product of all `bias` values (in biased motile states) and
the transition weights pass through every `transition_weights!`.

## Defining behaviors

Per-step logic without parameters or memory: a plain function, called as
`f(microbe, model)` right after translation.

```julia
bounce!(microbe, model) = position(microbe)[1] > 90 && microbe.vel[1] > 0 && (microbe.vel = -microbe.vel)
add_agent!(model; motility, behaviors = (bounce!,))
```

Stateless modulations of any hook: `Behavior`.

```@docs
Behavior
```

```julia
light(pos) = exp(-pos[2] / 50)
phobic = Behavior(bias = (microbe, model) -> 1 + 2 * light(position(microbe)))
add_agent!(model; motility, behaviors = (BrownBerg(), phobic))
```

Behaviors with memory: a struct plus methods for the hooks it uses.

```julia
@kwdef mutable struct Fatigue
    rate::Float64 = 0.01
    level::Float64 = 0.0
end
MicrobeAgents.affect!(f::Fatigue, microbe, model) = (f.level += f.rate * abmtimestep(model))
MicrobeAgents.speed_factor(f::Fatigue, microbe) = exp(-f.level)
```

## Multiple chemical fields

Any model property holding an `AbstractChemicalField` is a chemical field, and
`:chemicalfield` is the default one. Chemotaxis behaviors sense the default field
unless given a `field` keyword; use one behavior per field, each with its own
parameters and internal state:

```julia
model = StandardABM(Microbe{2}, space, dt;
    properties = Dict(:chemicalfield => attractant, :repellent => repellent))
add_agent!(model; motility,
    behaviors = (
        BrownBerg(gain = 660),
        BrownBerg(field = :repellent, gain = -400, receptor_binding_constant = 30),
    ))
```

Biases multiply across behaviors. For `BrownBerg` (`exp(-g S)`) this makes the
signals additive in log-rate; for `Xie` and `Brumley` (`1 + β s`) the combination is
a product. Naming a field the model does not have throws an `ArgumentError`
when the microbe is added.

Custom behaviors can read a named field through `concentration(microbe, model, key)`,
`gradient(...)`, `time_derivative(...)` and `diffusivity(...)`, and should call
`MicrobeAgents.check_field(model, key)` in `initialize!` so that a wrong name fails
when the microbe is added. Note that the default key was renamed from
`:chemoattractant` to `:chemicalfield`.

## Sharing measurements

`concentration`, `gradient`, `time_derivative` and `diffusivity`
are evaluated at most once per microbe per step *per field*, however many behaviors ask
for them. For other expensive quantities, define a sensor behavior placed
before the behaviors that use it; it stores the value in its own field, and
the others read it by name or through `findbehavior`.

Cached values refer to the position at the start of `affect_step!`. A behavior
that moves the microbe (e.g. a wall rule) should call
`MicrobeAgents.invalidate_field_cache!(model)` afterwards, so later behaviors
and `bias` hooks see fresh values.

```@docs
behaviors
findbehavior
```

## Chemotaxis models

```@docs
BrownBerg
Brumley
Celani
Xie
SonMenolascina
Chemokinesis
SpeedDependentTurnRate
SpeedDependentFlick
```

## Performance

`StandardABM(Microbe{D}, ...)` has a non-concrete agent type (Agents.jl warns,
and dispatch is dynamic). For a homogeneous population, declare the concrete
type; `N` must match the number of states of the motility.

```julia
bs = (chemotaxis = BrownBerg(),)
model = StandardABM(Microbe{2,2,typeof(bs)}, space, dt; properties)
add_agent!(model; motility = RunTumble([30.0], 0.67, Isotropic(2)), behaviors = bs)
```

For mixed populations, keep `Microbe{D}` and pass `warn = false` to silence
the warning.
