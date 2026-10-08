# Composable Behaviors Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace "one agent type per chemotaxis model" with a single `Microbe{D,N,B}` that carries a tuple of composable behaviors hooking into the stepping pipeline.

**Architecture:** Behaviors are plain values stored in `microbe.behaviors` (Tuple or NamedTuple). The stepping pipeline calls five hook functions (`initialize!`, `affect!`, `bias`, `speed_factor`, `transition_weights!`) on each behavior in tuple order, via compile-time-unrolled recursion over the tuple. Field quantities (`concentration`, `gradient`, ...) are memoized per microbe step in a model-level `FieldCache`. Existing chemotaxis models are rewritten as behavior structs and must reproduce pre-refactor trajectories bit-for-bit.

**Tech Stack:** Julia, Agents.jl 6/7, StaticArrays, StatsBase, Distributions, Documenter + Literate (docs).

**Spec:** `docs/superpowers/specs/2026-10-08-composable-behaviors-design.md`

## Global Constraints

- Work only inside the worktree `/home/riccardo/.julia/dev/MicrobeAgents/.claude/worktrees/behaviors-design` (branch `worktree-behaviors-design`). Other agents use the main checkout.
- Run Julia through the julia-mcp tools (`julia_eval` with `env_path` = worktree root), not `julia` on the command line (project directive). The session has no Revise: **after any change under `src/`, call `julia_restart` for that env before evaluating**.
- Agents.jl compat stays `"6, 7"`; no new package dependencies.
- Every exported function/type gets a short, direct docstring (no wordy LLM-style prose).
- `StandardABM(Microbe{D}, space, dt; properties)` keeps working unchanged; `add_agent!` only gains the `behaviors` keyword.
- Hook signatures, exactly: `initialize!(b, microbe, model)`, `affect!(b, microbe, model)`, `bias(b, microbe, model)`, `speed_factor(b, microbe)`, `transition_weights!(w, b, microbe, model)`. Plain `Function` behavior = `f(microbe, model)` as `affect!`.
- `bias` and `speed_factor` combine by product; others run in tuple order.
- Migrated chemotaxis models must reproduce the reference fixtures from Task 1 exactly (`==`).
- Commit after each task; messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Run git commands as plain single commands from the worktree root.

## Review Focus

1. **Single behavior passed without a tuple** (`behaviors = BrownBerg()`): should raise a clear `ArgumentError`, not an obscure MethodError deep in stepping. Test in Task 3.
2. **Behaviors shared across agents**: one tuple instance passed to many `add_agent!` calls must give each agent independent state, while plain functions/closures (which may capture large data) are *not* deep-copied. Test in Task 3.
3. **Noise-free sensing on a zero-radius microbe** (`chemotactic_precision = 0`, default `radius = 0`): must not produce `NaN`/throw (old code relied on `radius = 0.5` defaults; `0 * Inf = NaN`). Test in Tasks 5 and 6.
4. **Stale field cache** when a user reorders subroutines (e.g. `affect_step!` before `move_step!`) or calls `concentration` from `adata`: values must always correspond to the current position. Test in Task 4.
5. **Mixed populations** (agents with different behavior tuples, and with none) in one `StandardABM(Microbe{D}, ...)` with `container = Vector` and `Dict`. Test in Task 3.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/behaviors.jl` (new) | hook generics + defaults, `Behavior` wrapper, tuple combinators, `behaviors`, `findbehavior`, per-agent copying, `initialize_behaviors!` |
| `src/fields.jl` | + `FieldCache{D}`, reset/invalidate, memoized microbe-level wrappers |
| `src/microbes.jl` | `Microbe{D,N,B}`, partial constructor, `behaviors(::Microbe)`, `show` |
| `src/microbe_step.jl` | pipeline: cache handling, `affect_step!` over behaviors, aggregate `bias` in `switching_probability` |
| `src/motility.jl` | `update_motilestate!(microbe, model)` with scratch weights + `transition_weights!` hooks |
| `src/api.jl` | `speed(m)` applies speed factors; remove `state(::AbstractMicrobe)` |
| `src/model.jl` | default properties (`:field_cache`, no `:affect!`), `add_agent!` copies + initializes behaviors |
| `src/chemotaxis/sensing.jl` (new) | `CONV_NOISE`, `check_sensing_radius` |
| `src/chemotaxis/{brown-berg,brumley,celani,xie,son-menolascina}.jl` | rewritten as behaviors |
| `test/fixtures/generate_reference.jl`, `test/fixtures/reference.jl` (new) | pre-refactor reference outputs |
| `test/behaviors.jl`, `test/field-cache.jl`, `test/regression.jl` (new) | new tests |
| `test/{model-creation,model-stepping,chemotaxis,utils}.jl`, `test/runtests.jl` | migrated |
| `examples/Chemotaxis/*.jl`, `examples/Pathfinder/*.jl` | migrated |
| `docs/src/{introduction,api,behaviors}.md`, `docs/make.jl` | docs |

---

### Task 1: Reference fixtures from the current API

Must be done **before any change under `src/`**, on the unmodified code (commit `733ef48` + spec commits).

**Files:**
- Create: `test/fixtures/generate_reference.jl`
- Create (generated): `test/fixtures/reference.jl`

**Interfaces:**
- Produces: `test/fixtures/reference.jl` defining `const REFERENCE::Dict{String,<:NamedTuple}` keyed by `"BrownBerg"`, `"Brumley"`, `"Celani"`, `"Xie"`, `"SonMenolascina"`; each value has `Vector{Float64}` fields `xpos, ypos, speed, chemotaxis_bias, chemotaxis_state` (Xie also `state_m, state_z`), rows ordered as in the Agents.jl `adf` (time, then id), sampled at steps `0:10:300`.
- Produces: the scenario constants (`L`, `DT`, `NSTEPS`, `NAGENTS`, seed `1234`, fields `field/grad/ramp`) that Task 5–7's `test/regression.jl` must replicate exactly.

- [ ] **Step 1: Write the generator**

```julia
# Reference outputs of the chemotaxis models with the pre-behaviors API.
# Run on commit 733ef48 (before the composable-behaviors refactor), from the repo root:
#   using MicrobeAgents; include("test/fixtures/generate_reference.jl")
# It will NOT run on later versions; test/regression.jl replays the same scenarios
# with the behaviors API and compares against the generated reference.jl.
using MicrobeAgents, Random

const L = 1000.0
const DT = 0.1
const NSTEPS = 300
const NAGENTS = 3

field(m, model) = 1.0 + 0.01 * position(m)[1] + 0.001 * abmtime(model) * DT
grad(m::AbstractMicrobe{D}, model) where {D} = SVector{D}(i == 1 ? 0.01 : 0.0 for i in 1:D)
ramp(m, model) = 0.001
const CHEMO = GenericChemoattractant{2}(;
    concentration_field = field,
    concentration_gradient = grad,
    concentration_ramp = ramp,
)

xpos(m) = position(m)[1]
ypos(m) = position(m)[2]
chemotaxis_bias(m) = bias(m)
chemotaxis_state(m) = m.state
state_m(m) = m.state_m
state_z(m) = m.state_z

const RT = () -> RunTumble([30.0], 1.0, Isotropic(2); tumble_duration = 0.1)
const RRF = () -> RunReverseFlick([46.5], 0.45, [46.5], 0.45)

const SCENARIOS = [
    ("BrownBerg", BrownBerg{2}, RT, (rotational_diffusivity = 0.035, radius = 0.5), []),
    ("Brumley", Brumley{2}, RRF,
        (rotational_diffusivity = 0.035, radius = 0.5, chemotactic_precision = 6.0), []),
    ("Celani", Celani{2}, RT,
        (rotational_diffusivity = 0.26, radius = 0.5, chemotactic_precision = 6.0), []),
    ("Xie", Xie{2}, RRF,
        (rotational_diffusivity = 0.26, radius = 0.5, chemotactic_precision = 6.0),
        [state_m, state_z]),
    ("SonMenolascina", SonMenolascina{2},
        () -> RunReverseFlick([30.0], 0.5, [30.0], 0.5),
        (rotational_diffusivity = 0.035, radius = 0.5), []),
]

function run_scenario(T, motility, kw, extras)
    space = ContinuousSpace((L, L); periodic = true)
    model = StandardABM(T, space, DT;
        properties = Dict(:chemoattractant => CHEMO),
        rng = Xoshiro(1234), container = Vector,
    )
    for _ in 1:NAGENTS
        add_agent!(model; motility = motility(), kw...)
    end
    adata = [xpos, ypos, speed, chemotaxis_bias, chemotaxis_state, extras...]
    adf, = run!(model, NSTEPS; adata, when = 0:10:NSTEPS)
    cols = Symbol.(names(adf))
    keep = filter(c -> c ∉ (:time, :id), cols)
    NamedTuple{Tuple(keep)}(Tuple(Vector{Float64}(adf[!, c]) for c in keep))
end

open(joinpath(@__DIR__, "reference.jl"), "w") do io
    println(io, "# Generated by generate_reference.jl on commit 733ef48. Do not edit.")
    println(io, "const REFERENCE = Dict(")
    for (name, T, motility, kw, extras) in SCENARIOS
        cols = run_scenario(T, motility, kw, extras)
        println(io, "    ", repr(name), " => (")
        for (k, v) in pairs(cols)
            println(io, "        ", k, " = ", repr(v), ",")
        end
        println(io, "    ),")
    end
    println(io, ")")
end
```

- [ ] **Step 2: Run it on the unmodified code**

Confirm `git diff 733ef48 -- src` is empty first. Then julia-mcp, `env_path` = worktree root:
`using MicrobeAgents; include("test/fixtures/generate_reference.jl")`
Expected: `test/fixtures/reference.jl` exists (~50–80 KB).

- [ ] **Step 3: Sanity-check the fixture**

julia-mcp: `include("test/fixtures/reference.jl"); for (k, v) in REFERENCE; println(k, " ", keys(v), " ", length(v.xpos), " ", any(isnan, v.chemotaxis_bias)); end`
Expected: 5 lines, each with 93 rows (31 samples × 3 agents), `false` for NaN; Xie has `state_m`, `state_z`.

- [ ] **Step 4: Commit**

```bash
git add test/fixtures/generate_reference.jl test/fixtures/reference.jl
git commit -m "test: add pre-refactor reference outputs for chemotaxis models"
```

---

### Task 2: Behavior hooks, `Behavior` wrapper, combinators

**Files:**
- Create: `src/behaviors.jl`
- Modify: `src/MicrobeAgents.jl` (include after `fields.jl`, before `microbes.jl`)
- Modify: `src/microbes.jl` (remove `bias` from its export list only — `bias` is now exported from `behaviors.jl`; keep the old 1-arg `bias(microbe::AbstractMicrobe) = 1.0` method until Task 3 since `switching_probability` still calls it)
- Test: `test/behaviors.jl`

**Interfaces:**
- Produces (exported): `initialize!`, `affect!`, `bias`, `speed_factor`, `transition_weights!`, `Behavior`, `behaviors`, `findbehavior`.
- Produces (internal): `_foreach(f, t::Tuple)`, `_prod(f, t::Tuple)::Float64`, `copybehavior(b)`, `copybehaviors(bs)` (Tuple/NamedTuple → same container type; anything else → `ArgumentError`), `initialize_behaviors!(microbe, model)`, `behaviors(m::AbstractMicrobe) = ()`.
- Aggregate: `bias(microbe::AbstractMicrobe, model)` = product of `bias(b, microbe, model)`.

- [ ] **Step 1: Write the failing tests**

```julia
using MicrobeAgents, Test
using MicrobeAgents: _foreach, _prod, copybehavior, copybehaviors

# types and hook methods must be defined at top level, not inside @testset
struct Dummy end
mutable struct Counter
    n::Int
end

@testset "Behaviors" begin
    @testset "Default hooks are neutral" begin
        @test isnothing(initialize!(Dummy(), nothing, nothing))
        @test isnothing(affect!(Dummy(), nothing, nothing))
        @test bias(Dummy(), nothing, nothing) === 1.0
        @test speed_factor(Dummy(), nothing) === 1.0
        w = [0.5, 0.5]
        @test isnothing(transition_weights!(w, Dummy(), nothing, nothing))
        @test w == [0.5, 0.5]
    end

    @testset "Plain functions are affect!" begin
        hits = Ref(0)
        f(microbe, model) = (hits[] += 1; :ignored)
        @test isnothing(affect!(f, nothing, nothing))
        @test hits[] == 1
        @test bias(f, nothing, nothing) === 1.0
    end

    @testset "Behavior wrapper" begin
        b = Behavior(
            bias = (m, model) -> 2.0,
            speed_factor = m -> 3.0,
            transition_weights! = (w, m, model) -> (w[1] = 0.0),
        )
        @test bias(b, nothing, nothing) === 2.0
        @test speed_factor(b, nothing) === 3.0
        w = [0.5, 0.5]
        transition_weights!(w, b, nothing, nothing)
        @test w == [0.0, 0.5]
        @test isnothing(affect!(b, nothing, nothing))
        hits = Ref(0)
        b2 = Behavior(affect! = (m, model) -> (hits[] += 1))
        affect!(b2, nothing, nothing)
        @test hits[] == 1
        @test bias(b2, nothing, nothing) === 1.0
        @test_throws MethodError Behavior(unknown = identity)
    end

    @testset "Combinators" begin
        order = Int[]
        _foreach(x -> push!(order, x), (1, 2, 3))
        @test order == [1, 2, 3]
        @test isnothing(_foreach(identity, ()))
        @test _prod(identity, ()) === 1.0
        @test _prod(identity, (2.0, 3.0)) === 6.0
        @test (@inferred _prod(x -> x isa Int ? 2.0 : 0.5, (1, "a", 3))) === 2.0
    end

    @testset "Copying" begin
        f(m, model) = nothing
        c = Counter(0)
        bs = (c, f)
        cs = copybehaviors(bs)
        @test cs[1] !== c && cs[1].n == 0
        @test cs[2] === f
        nt = (a = Counter(1), b = Behavior(bias = (m, model) -> 2.0))
        cnt = copybehaviors(nt)
        @test cnt isa NamedTuple{(:a, :b)}
        @test cnt.a !== nt.a && cnt.b === nt.b
        @test_throws ArgumentError copybehaviors(Counter(0))
    end
end
```

- [ ] **Step 2: Run tests to verify they fail**

julia-mcp (env = worktree): `using MicrobeAgents, Test; include("test/behaviors.jl")`
Expected: errors — `initialize!`/`Behavior` not defined.

- [ ] **Step 3: Implement `src/behaviors.jl`**

```julia
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
```

In `src/MicrobeAgents.jl` add `include("behaviors.jl")` right after `include("fields.jl")`.

- [ ] **Step 4: Run tests to verify they pass**

`julia_restart`, then: `using MicrobeAgents, Test; include("test/behaviors.jl")`
Expected: all pass. (The rest of the suite is still on the old API; old chemotaxis files still define `bias(microbe::BrownBerg)` etc. as 1-arg methods of the same generic, which is harmless for now.)

- [ ] **Step 5: Commit**

```bash
git add src/behaviors.jl src/MicrobeAgents.jl src/microbes.jl test/behaviors.jl
git commit -m "feat: add behavior hooks, Behavior wrapper and tuple combinators"
```

---

### Task 3: `Microbe{D,N,B}`, `add_agent!`, stepping pipeline

This task removes the old agent-type chemotaxis models from the build (re-added as behaviors in Tasks 5–7).

**Files:**
- Modify: `src/microbes.jl` (struct, constructor, accessor, show; drop `chemotaxis!`)
- Modify: `src/microbe_step.jl` (`affect_step!`, `switching_probability`)
- Modify: `src/motility.jl:207-217` (`update_motilestate!(microbe, model)`)
- Modify: `src/api.jl` (`speed`, remove `state(m::AbstractMicrobe)` + its docstring and `state` from api.jl's export line — `state` stays exported from `motility.jl` for `Motility`)
- Modify: `src/model.jl` (`add_agent!`, default properties, docstring)
- Modify: `src/MicrobeAgents.jl` (comment out the five `include("chemotaxis/...")` lines; update `AbstractMicrobe` docstring field list: `id, pos, vel, speed, motility, rotational_diffusivity, radius`, plus optional `behaviors(m)` method)
- Modify: `test/runtests.jl`, `test/model-creation.jl`, `test/model-stepping.jl`, `test/utils.jl`
- Test: `test/behaviors.jl` (extend)

**Interfaces:**
- Consumes: everything from Task 2.
- Produces: `Microbe{D,N,B}` with fields `id, pos, vel, speed, motility, rotational_diffusivity, radius, behaviors`; `Microbe{D,N}(; behaviors=(), kw...)`; `behaviors(m::Microbe)`; `add_agent!(...; motility, behaviors=(), ...)` returning the microbe; `speed(m) = m.speed * _prod(b -> speed_factor(b, m), values(behaviors(m)))`; `affect_step!` loops `affect!`; `switching_probability` uses `bias(microbe, model)`; `update_motilestate!(microbe, model)` applies `transition_weights!` on a scratch `ProbabilityWeights`.
- Default model properties: `:chemoattractant`, `:timestep` (Task 4 adds `:field_cache`).

- [ ] **Step 1: Write failing tests**

Add at the top level of `test/behaviors.jl` (next to `Dummy`/`Counter`):
```julia
using Random
mutable struct Ticker
    n::Int
end
MicrobeAgents.affect!(t::Ticker, microbe, model) = (t.n += 1; nothing)
MicrobeAgents.initialize!(t::Ticker, microbe, model) = (t.n = 100; nothing)
```
Then append inside the top-level `@testset "Behaviors"`:

```julia
    @testset "Microbe with behaviors" begin
        for D in 1:3, container in (Vector, Dict)
            space = ContinuousSpace(fill(100.0, SVector{D}))
            model = StandardABM(Microbe{D}, space, 0.1; container, rng = Xoshiro(1))
            motility = RunTumble([30.0], 1.0, Isotropic(D))
            shared = (t = Ticker(0),)
            add_agent!(model; motility, behaviors = shared)
            add_agent!(model; motility, behaviors = shared)
            add_agent!(model; motility) # no behaviors: mixed population
            @test model[1] isa Microbe{D,2}
            @test behaviors(model[3]) === ()
            @test model[1].behaviors.t !== model[2].behaviors.t # independent copies
            @test shared.t.n == 0 # original untouched
            @test model[1].behaviors.t.n == 100 # initialize! ran
            run!(model, 3)
            @test model[1].behaviors.t.n == model[2].behaviors.t.n == 103
            @test findbehavior(model[1], Ticker) === model[1].behaviors.t
            @test findbehavior(model[3], Ticker) === nothing
        end
    end

    @testset "Single behavior without tuple" begin
        space = ContinuousSpace((10.0, 10.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        motility = RunTumble([30.0], 1.0, Isotropic(2))
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = Behavior(bias = (m, model) -> 2.0))
    end

    @testset "Hooks in the pipeline" begin
        space = ContinuousSpace((1000.0, 1000.0))
        dt = 0.1
        # bias: product over behaviors, used in switching_probability
        model = StandardABM(Microbe{2}, space, dt)
        motility = RunTumble([30.0], 1.0, Isotropic(2))
        b2 = Behavior(bias = (m, model) -> 2.0)
        b3 = Behavior(bias = (m, model) -> 3.0)
        add_agent!(model; motility, behaviors = (b2, b3))
        @test bias(model[1], model) == 6.0
        @test switching_probability(model[1], model) ≈ 6.0 * dt / 1.0
        # speed_factor: product, applied on top of the sampled speed
        add_agent!(model; motility, speed = 10.0,
            behaviors = (Behavior(speed_factor = m -> 2.0), Behavior(speed_factor = m -> 1.5)))
        @test model[2].speed == 10.0
        @test speed(model[2]) == 30.0
        @test velocity(model[2]) ≈ direction(model[2]) .* 30.0
        # plain function as affect!
        hits = Ref(0)
        counter(m, model) = (hits[] += 1)
        add_agent!(model; motility, behaviors = (counter,))
        run!(model, 2)
        @test hits[] == 2
    end

    @testset "transition_weights! uses a scratch copy" begin
        space = ContinuousSpace((1000.0, 1000.0))
        model = StandardABM(Microbe{2}, space, 1.0; rng = Xoshiro(3))
        # RunReverseFlick: state 3 (backward run) always goes to 4 (flick);
        # the hook redirects it to state 2 (reverse) instead
        motility = RunReverseFlick([30.0], 0.0, [30.0], 0.0)
        redirect = Behavior(transition_weights! = (w, m, model) ->
            (state(motilepattern(m)) == 3 && (w[2] = 1.0; w[4] = 0.0)))
        add_agent!(model; motility, behaviors = (redirect,))
        m = model[1]
        m.motility.current_state = 3
        MicrobeAgents.update_motilestate!(m, model)
        @test state(motilepattern(m)) == 2
        @test transition_weights(motilepattern(m), 3) == [0.0, 0.0, 0.0, 1.0] # unchanged
    end
```

- [ ] **Step 2: Run to verify failure**

julia-mcp: `using MicrobeAgents, Test; include("test/behaviors.jl")`
Expected: failures/errors (`behaviors` keyword not accepted / no field `behaviors`).

- [ ] **Step 3: Rewrite `src/microbes.jl`**

```julia
export Microbe

"""
    Microbe{D,N,B} <: AbstractMicrobe{D,N}
Microbe type of MicrobeAgents.jl.
`D` is the space dimensionality, `N` the number of motile states,
`B` the type of the `behaviors` container (inferred, never written by hand).

Fields:
- `id::Int` identifier used internally
- `pos::SVector{D,Float64}` position
- `vel::SVector{D,Float64}` unit direction of motion
- `speed::Float64` sampled speed of the current motile state (see `speed`)
- `motility::Motility{N}` motility pattern
- `rotational_diffusivity::Float64 = 0.0` rotational diffusion coefficient
- `radius::Float64 = 0.0` equivalent spherical radius
- `behaviors::B = ()` `Tuple` or `NamedTuple` of behaviors (see `Behavior`)
"""
@agent struct Microbe{D,N,B}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N}
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    behaviors::B = ()
end
Microbe{D,N}(; behaviors = (), kwargs...) where {D,N} =
    Microbe{D,N,typeof(behaviors)}(; behaviors, kwargs...)

behaviors(m::Microbe) = m.behaviors

r2dig(x) = round(x, digits=2)
_label(b) = string(nameof(typeof(b)))
_label(f::Function) = string(nameof(f))
_labels(bs::Tuple) = map(_label, bs)
_labels(bs::NamedTuple) = map((k, b) -> "$k = $(_label(b))", keys(bs), values(bs))
function Base.show(io::IO, ::MIME"text/plain", m::AbstractMicrobe{D,N}) where {D,N}
    println(io, "$(nameof(typeof(m))){$D} with $(N)-state motility pattern")
    println(io, "position (μm): $(r2dig.(position(m))); velocity (μm/s): $(r2dig.(velocity(m)))")
    bs = behaviors(m)
    print(io, "behaviors: ", isempty(bs) ? "none" : join(_labels(bs), ", "))
end
```

Verify in julia-mcp after restart that `Microbe{2,2}(; id=1, pos=SVector(0.0,0.0), vel=SVector(1.0,0.0), speed=1.0, motility=RunTumble([1.0],1.0,Isotropic(2)))` returns a `Microbe{2,2,Tuple{}}`. If `@agent` does not generate the `Microbe{D,N,B}(; ...)` keyword constructor, add it explicitly:
```julia
Microbe{D,N,B}(; id, pos, vel, speed, motility, rotational_diffusivity = 0.0,
    radius = 0.0, behaviors) where {D,N,B} =
    Microbe{D,N,B}(id, pos, vel, speed, motility, rotational_diffusivity, radius, behaviors)
```

Also delete from `src/microbes.jl` the old `chemotaxis!` fallback and the 1-arg `bias(microbe::AbstractMicrobe) = 1.0`, with their docstrings.

- [ ] **Step 4: Update `src/api.jl`**

Replace the `speed` method and docstring:
```julia
"""
    speed(m::AbstractMicrobe)
Return the speed of the microbe: the speed sampled from its current motile
state, times the `speed_factor` of all its behaviors.
"""
speed(m::AbstractMicrobe) = m.speed * _prod(b -> speed_factor(b, m), values(behaviors(m)))
```
Delete the `state(m::AbstractMicrobe)` method with its docstring, and remove `state` from the `export` line at the top of `api.jl`.

- [ ] **Step 5: Update `src/microbe_step.jl`**

Replace `affect_step!` and the `β` line in `switching_probability`:
```julia
"""
    affect_step!(microbe::AbstractMicrobe, model::ABM)
Subroutine calling `affect!` for each behavior of the microbe, in order.
"""
function affect_step!(microbe::AbstractMicrobe, model::ABM)
    _foreach(b -> affect!(b, microbe, model), values(behaviors(microbe)))
end
```
```julia
    β = biased(M) ? bias(microbe, model) : 1.0
```
Update the `switching_probability` docstring bullet to: "β is the total `bias(microbe, model)` of the microbe behaviors (product of their `bias` hooks); β=1 is unbiased." Update the `microbe_step!` docstring step 3 to "Call the `affect!` hook of each behavior".

- [ ] **Step 6: Update `src/motility.jl` `update_motilestate!(microbe, model)`**

```julia
"""
    update_motilestate!(microbe, model)
Update the motile state of `microbe` by randomly sampling the next state
according to the transition weights, after the `transition_weights!`
hooks of its behaviors have been applied to a copy of them.
"""
function update_motilestate!(microbe::AbstractMicrobe, model::AgentBasedModel)
    motility = motilepattern(microbe)
    w0 = transition_weights(motility, state(motility))
    w = ProbabilityWeights(copy(w0.values), sum(w0))
    _foreach(b -> transition_weights!(w, b, microbe, model), values(behaviors(microbe)))
    j = sample(abmrng(model), eachindex(w), w)
    update_motilestate!(motility, j)
end
```
(Keep `update_motilestate!(motility::Motility, model)` and `update_motilestate!(motility::Motility, j::Int)` unchanged.)

- [ ] **Step 7: Update `src/model.jl`**

In the main `add_agent!` (the `Type{<:AbstractMicrobe{D,N}}` method):
```julia
    id = Agents.nextid(model) # not public API!
    if !isempty(properties)
        microbe = A(id, pos, properties...)
    else
        kw = haskey(kwproperties, :behaviors) ?
            merge(values(kwproperties), (behaviors = copybehaviors(kwproperties[:behaviors]),)) :
            kwproperties
        microbe = A(; id, pos, vel = zero(SVector{D}), speed = 0.0, kw...)
        microbe.vel = isnothing(vel) ? random_velocity(model) : vel
        microbe.speed = isnothing(speed) ? random_speed(microbe, model) : speed
    end
    # initialize before placement: a throwing initialize! leaves the model untouched
    initialize_behaviors!(microbe, model)
    Agents.add_agent_own_pos!(microbe, model) # not public API!
    return microbe
```
Extend its docstring: "`behaviors` (a `Tuple` or `NamedTuple`) is copied for each agent (functions and `Behavior`s are shared), and `initialize!` is called on each behavior before the microbe is placed in the model."
Change `make_default_abm_properties(D)` to `Dict(:chemoattractant => GenericChemoattractant{D}())` and the `StandardABM` docstring's default-properties block accordingly (remove `:affect! => chemotaxis!`).

- [ ] **Step 8: Drop old chemotaxis from build and migrate old tests**

- `src/MicrobeAgents.jl`: comment out the five `include("chemotaxis/...")` lines with `# re-enabled as behaviors in later tasks`.
- `test/runtests.jl`: add `include("behaviors.jl")` after `utils.jl`; comment out `include("chemotaxis.jl")`.
- `test/model-creation.jl`: property set becomes `Set((:timestep, :chemoattractant))`; remove `@test state(model[1]) == 0.0`; delete the whole `"Chemotactic Microbe types"` testset (replaced in Tasks 5–7).
- `test/model-stepping.jl`: in both places, replace the custom `affect!` property block by
```julia
            # custom per-step logic as a behavior
            mutable struct Countdown
                s::Float64
            end
            MicrobeAgents.affect!(c::Countdown, microbe::Microbe{D}, model) where {D} = (c.s -= D)
            model = StandardABM(Microbe{D}, space, dt; container)
            motility = RunTumble(
                run_duration=1.0, run_speed=[30.0],
                angle=Isotropic(D), tumble_duration=0.0
            )
            add_agent!(model; motility, behaviors = (Countdown(0.0),))
            run!(model, 1)
            @test model[1].behaviors[1].s == -D
```
  Move the `mutable struct Countdown` and the `MicrobeAgents.affect!(c::Countdown, ...)` method out of the loop to the top level of the file (below the `using` lines): type and qualified method definitions are not allowed in local scope.
- `test/utils.jl`: `MicrobeTypes = [Microbe{D}]` for now (Task 5 restores variety via behaviors).

- [ ] **Step 9: Run tests**

`julia_restart`, then `using Pkg; Pkg.test()`
Expected: all included testsets pass.

- [ ] **Step 10: Commit**

```bash
git add -A src test
git commit -m "feat!: Microbe carries composable behaviors; pipeline calls behavior hooks"
```

---

### Task 4: Memoized field quantities

**Files:**
- Modify: `src/fields.jl`
- Modify: `src/microbe_step.jl` (cache reset/invalidate)
- Modify: `src/model.jl` (`:field_cache` default property + docstring)
- Modify: `test/model-creation.jl` (property set)
- Test: `test/field-cache.jl` (new), add to `test/runtests.jl` after `behaviors.jl`

**Interfaces:**
- Produces (internal): `FieldCache{D}`, `field_cache(model)`, `reset_field_cache!(model, microbe)`, `invalidate_field_cache!(model)`.
- Behavior: microbe-level `concentration/gradient/time_derivative/chemoattractant_diffusivity(microbe, model)` return memoized values when called for the microbe currently in its `affect!`/reorient phase; otherwise compute fresh.

- [ ] **Step 1: Write failing tests (`test/field-cache.jl`)**

```julia
using MicrobeAgents, Test, Random

@testset "Field cache" begin
    calls = Ref(0)
    cfield(m, model) = (calls[] += 1; position(m)[1])
    chemo = GenericChemoattractant{2}(; concentration_field = cfield)
    space = ContinuousSpace((100.0, 100.0))
    model = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemoattractant => chemo))
    @test :field_cache in keys(abmproperties(model))
    reader(m, model) = concentration(m, model)
    motility = RunTumble([1.0], Inf, Isotropic(2))
    add_agent!(SVector(10.0, 50.0), model; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader, reader, reader))
    calls[] = 0
    run!(model, 1)
    @test calls[] == 1 # three readers, one evaluation

    # outside a step: always fresh
    calls[] = 0
    @test concentration(model[1], model) == 11.0
    @test concentration(model[1], model) == 11.0
    @test calls[] == 2

    # reordered subroutines: affect before move must not leak stale values into reorient
    seen = Float64[]
    probe = Behavior(bias = (m, model) -> (push!(seen, concentration(m, model)); 1.0))
    step_affect_first!(m, model) = (affect_step!(m, model); move_step!(m, model); reorient_step!(m, model))
    model2 = StandardABM(Microbe{2}, space, 1.0;
        properties = Dict(:chemoattractant => chemo), agent_step! = step_affect_first!)
    add_agent!(SVector(10.0, 50.0), model2; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader, probe))
    run!(model2, 1)
    @test seen == [11.0] # post-move position, not the cached pre-move 10.0

    # type stability
    m = model[1]
    @test (@inferred concentration(m, model)) isa Float64
    @test (@inferred gradient(m, model)) isa SVector{2,Float64}
end
```

- [ ] **Step 2: Run to verify failure**

`julia_restart`; `using MicrobeAgents, Test; include("test/field-cache.jl")`
Expected: FAIL (`:field_cache` missing, `calls[] == 3`).

- [ ] **Step 3: Implement cache in `src/fields.jl`**

Replace the four microbe-level wrappers with:
```julia
# per-step memoization of field quantities for the microbe being stepped;
# `id == 0` means no microbe is being stepped
mutable struct FieldCache{D}
    id::Int
    concentration::Union{Nothing,Float64}
    gradient::Union{Nothing,SVector{D,Float64}}
    time_derivative::Union{Nothing,Float64}
    chemoattractant_diffusivity::Union{Nothing,Float64}
end
FieldCache{D}() where {D} = FieldCache{D}(0, nothing, nothing, nothing, nothing)

field_cache(model::ABM) = abmproperties(model).field_cache

function reset_field_cache!(model::ABM, microbe::AbstractMicrobe)
    c = field_cache(model)
    c.id = microbe.id
    c.concentration = nothing
    c.gradient = nothing
    c.time_derivative = nothing
    c.chemoattractant_diffusivity = nothing
    return nothing
end
invalidate_field_cache!(model::ABM) = (field_cache(model).id = 0; nothing)

@inline function _cached(compute::F, name::Symbol, microbe, model) where {F}
    cache = field_cache(model)
    cache.id == microbe.id || return compute()
    v = getfield(cache, name)
    v === nothing || return v
    v = compute()
    setfield!(cache, name, v)
    return v
end

function concentration(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:concentration, microbe, model) do
        concentration(chemoattractant(model))(microbe, model)::Float64
    end
end
function gradient(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:gradient, microbe, model) do
        gradient(chemoattractant(model))(microbe, model)::SVector{D,Float64}
    end
end
function time_derivative(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:time_derivative, microbe, model) do
        time_derivative(chemoattractant(model))(microbe, model)::Float64
    end
end
function chemoattractant_diffusivity(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:chemoattractant_diffusivity, microbe, model) do
        chemoattractant_diffusivity(chemoattractant(model))(microbe, model)::Float64
    end
end
```
If `@inferred` fails in Step 4, replace `getfield(cache, name)`/`setfield!` with four explicit per-field branches (no Symbol argument).

- [ ] **Step 4: Wire into the pipeline (`src/microbe_step.jl`)**

- First line of `move_step!` and of `move_step_pathfinder!`: `invalidate_field_cache!(model)`.
- First line of `affect_step!`: `reset_field_cache!(model, microbe)`; docstring: "Resets the per-step cache of field quantities, then calls `affect!` for each behavior of the microbe, in order."
- Last line of `reorient_step!`: `invalidate_field_cache!(model)` followed by `return nothing`.

In `src/model.jl`: `make_default_abm_properties(D) = Dict(:chemoattractant => GenericChemoattractant{D}(), :field_cache => FieldCache{D}())`, and list `:field_cache => FieldCache{D}()` (internal, per-step cache of field quantities) in the `StandardABM` docstring.
In `test/model-creation.jl`: `Set((:timestep, :chemoattractant, :field_cache))`.
In `test/runtests.jl`: `include("field-cache.jl")` after `include("behaviors.jl")`.

- [ ] **Step 5: Run tests**

`julia_restart`; `using Pkg; Pkg.test()`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -A src test
git commit -m "feat: memoize field quantities once per microbe step"
```

---

### Task 5: BrownBerg and Brumley as behaviors + regression harness

**Files:**
- Create: `src/chemotaxis/sensing.jl`
- Rewrite: `src/chemotaxis/brown-berg.jl`, `src/chemotaxis/brumley.jl`
- Modify: `src/MicrobeAgents.jl` (move `CONV_NOISE` into `sensing.jl`; include `chemotaxis/sensing.jl`, `brown-berg.jl`, `brumley.jl`)
- Create: `test/regression.jl`; add `include("regression.jl")` to `test/runtests.jl` (after `field-cache.jl`)
- Rewrite: `test/chemotaxis.jl` (BrownBerg, Brumley sections now; re-enable its include in `runtests.jl`)
- Modify: `test/utils.jl` (distance test over behaviors instead of types)

**Interfaces:**
- Consumes: hooks (Task 2), `Microbe`/`add_agent!` (Task 3), cached `concentration` etc. (Task 4).
- Produces: `BrownBerg(; gain=660.0, receptor_binding_constant=100.0, memory=1.0)` mutable, field `state`; `Brumley(; memory=1.3, gain_receptor=50.0, gain=50.0, chemotactic_precision=6.0)` mutable, field `state`; internal `check_sensing_radius(behavior, Π, microbe)`; `CONV_NOISE`.
- `test/regression.jl` exposes `run_scenario(motility, behaviors, kw, extras)` used by Tasks 6–7.

- [ ] **Step 1: Write the regression harness (`test/regression.jl`)**

```julia
using MicrobeAgents, Test, Random
include(joinpath(@__DIR__, "fixtures", "reference.jl"))

# same scenarios as fixtures/generate_reference.jl, with the behaviors API
const REF_L = 1000.0
const REF_DT = 0.1
const REF_NSTEPS = 300
const REF_NAGENTS = 3
ref_field(m, model) = 1.0 + 0.01 * position(m)[1] + 0.001 * abmtime(model) * REF_DT
ref_grad(m::AbstractMicrobe{D}, model) where {D} = SVector{D}(i == 1 ? 0.01 : 0.0 for i in 1:D)
ref_ramp(m, model) = 0.001
const REF_CHEMO = GenericChemoattractant{2}(;
    concentration_field = ref_field,
    concentration_gradient = ref_grad,
    concentration_ramp = ref_ramp,
)
ref_rt() = RunTumble([30.0], 1.0, Isotropic(2); tumble_duration = 0.1)
ref_rrf() = RunReverseFlick([46.5], 0.45, [46.5], 0.45)

function run_scenario(motility, behaviors, kw, extras)
    space = ContinuousSpace((REF_L, REF_L); periodic = true)
    model = StandardABM(Microbe{2}, space, REF_DT;
        properties = Dict(:chemoattractant => REF_CHEMO),
        rng = Xoshiro(1234), container = Vector,
    )
    for _ in 1:REF_NAGENTS
        add_agent!(model; motility = motility(), behaviors, kw...)
    end
    xpos(m) = position(m)[1]
    ypos(m) = position(m)[2]
    chemotaxis_bias(m) = bias(m.behaviors.chemotaxis, m, model)
    chemotaxis_state(m) = m.behaviors.chemotaxis.state
    adata = [xpos, ypos, speed, chemotaxis_bias, chemotaxis_state, extras...]
    adf, = run!(model, REF_NSTEPS; adata, when = 0:10:REF_NSTEPS)
    return adf
end

function test_against_reference(name, adf)
    ref = REFERENCE[name]
    @testset "$name reproduces reference" begin
        for k in keys(ref)
            @test Vector{Float64}(adf[!, k]) == ref[k]
        end
    end
end

@testset "Regression vs pre-behaviors API" begin
    test_against_reference("BrownBerg", run_scenario(ref_rt,
        (chemotaxis = BrownBerg(),),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
    test_against_reference("Brumley", run_scenario(ref_rrf,
        (chemotaxis = Brumley(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
end
```

- [ ] **Step 2: Run to verify failure**

`julia_restart`; `using MicrobeAgents, Test; include("test/regression.jl")`
Expected: `UndefVarError: BrownBerg`.

- [ ] **Step 3: Create `src/chemotaxis/sensing.jl`**

```julia
"""
Conversion factor (1/√(number of molecules) --> 1/√(moles)) used
in the evaluation of chemotactic sensing noise.
"""
global const CONV_NOISE::Float64 = 0.04075

# Berg-Purcell noise scales as 1/√radius, so noisy sensing needs a finite radius
function check_sensing_radius(behavior, Π, microbe)
    if Π > 0 && iszero(radius(microbe))
        throw(ArgumentError(
            "$(nameof(typeof(behavior))) with chemotactic_precision > 0 requires " *
            "a microbe radius > 0 (e.g. `add_agent!(model; radius = 0.5, ...)`)"
        ))
    end
    return nothing
end
```
In `src/MicrobeAgents.jl`, remove the `CONV_NOISE` docstring+definition and replace the chemotaxis block with:
```julia
# implementations of chemotactic models as behaviors
include("chemotaxis/sensing.jl")
include("chemotaxis/brown-berg.jl")
include("chemotaxis/brumley.jl")
# include("chemotaxis/celani.jl") # re-enabled in later tasks
# include("chemotaxis/xie.jl")
# include("chemotaxis/son-menolascina.jl")
```

- [ ] **Step 4: Rewrite `src/chemotaxis/brown-berg.jl`**

```julia
export BrownBerg

"""
    BrownBerg(; gain=660, receptor_binding_constant=100, memory=1)
Chemotaxis behavior from 'Brown and Berg (1974) PNAS'.

Parameters:
- `gain = 660` s
- `receptor_binding_constant = 100` μM
- `memory = 1` s
Internal state: `state`, the weighted dPb/dt of the paper.
"""
@kwdef mutable struct BrownBerg
    gain::Float64 = 660.0
    receptor_binding_constant::Float64 = 100.0
    memory::Float64 = 1.0
    state::Float64 = 0.0
end

function affect!(b::BrownBerg, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    τₘ = b.memory
    β = exp(-Δt / τₘ) # memory loss factor
    KD = b.receptor_binding_constant
    S = b.state # weighted dPb/dt at previous step
    vel = velocity(microbe)
    u = concentration(microbe, model)
    ∇u = gradient(microbe, model)
    ∂ₜu = time_derivative(microbe, model)
    du_dt = dot(vel, ∇u) + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    b.state = (1 - β) * M + S * β # new weighted dPb/dt
    return nothing
end

bias(b::BrownBerg, microbe, model) = exp(-b.gain * b.state)
```

- [ ] **Step 5: Rewrite `src/chemotaxis/brumley.jl`**

```julia
export Brumley

"""
    Brumley(; memory=1.3, gain_receptor=50, gain=50, chemotactic_precision=6)
Chemotaxis behavior from 'Brumley et al. (2019) PNAS', with gaussian sensing
noise in the gradient measurement. Requires a microbe `radius > 0` when
`chemotactic_precision > 0`.

Parameters:
- `memory = 1.3` s → 'τₘ'
- `gain_receptor = 50` μM⁻¹ → 'κ'
- `gain = 50` → 'Γ'
- `chemotactic_precision = 6` → 'Π'
Internal state: `state` → 'S'.
"""
@kwdef mutable struct Brumley
    memory::Float64 = 1.3
    gain_receptor::Float64 = 50.0
    gain::Float64 = 50.0
    chemotactic_precision::Float64 = 6.0
    state::Float64 = 0.0
end

initialize!(b::Brumley, microbe, model) =
    check_sensing_radius(b, b.chemotactic_precision, microbe)

function affect!(b::Brumley, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = chemoattractant_diffusivity(microbe, model)
    τₘ = b.memory
    α = exp(-Δt / τₘ) # memory persistence factor
    a = radius(microbe)
    Π = b.chemotactic_precision
    κ = b.gain_receptor
    vel = velocity(microbe)
    u = concentration(microbe, model)
    ∇u = gradient(microbe, model)
    ∂ₜu = time_derivative(microbe, model)
    # gradient measurement
    μ = dot(vel, ∇u) + ∂ₜu # mean
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * u / (π * a * Dc * Δt^3)) # noise
    M = rand(abmrng(model), Normal(μ, σ)) # measurement
    # update internal state
    S = b.state
    b.state = α * S + (1 - α) * κ * τₘ * M
    return nothing
end

bias(b::Brumley, microbe, model) = (1 + exp(-b.gain * b.state)) / 2
```

- [ ] **Step 6: Rewrite BrownBerg/Brumley sections of `test/chemotaxis.jl`**

Keep the helper field functions at the top. Convert each block mechanically:
`StandardABM(BrownBerg{2,2}, space, dt; properties)` → `StandardABM(Microbe{2}, space, dt; properties)`;
`add_agent!(model; motility, gain=600, receptor_binding_constant=100, memory=1)` → `add_agent!(model; motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=100, memory=1),))` (same for the `pos`/`vel` variants and Brumley with `Brumley(gain=..., memory=..., chemotactic_precision=0)`);
`bias(model[i])` → `bias(model[i], model)`;
`adata=[bias]` → define `tumblebias(m) = bias(m, model)` before `run!` and use `adata=[tumblebias]`, `Analysis.adf_to_vectors(adf, :tumblebias)`.
Comment out the Celani and SonMenolascina testsets (restored in Tasks 6–7). Add:
```julia
    @testset "Noise-free sensing with zero radius" begin
        space = ContinuousSpace((100.0, 100.0))
        chemo = GenericChemoattractant{2}(; concentration_field = constant_background_concentration)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:chemoattractant => chemo))
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = (Brumley(chemotactic_precision = 0),))
        run!(model, 5)
        @test bias(model[1], model) == 1
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = (Brumley(chemotactic_precision = 6),))
    end
```
Re-enable `include("chemotaxis.jl")` in `test/runtests.jl`.
In `test/utils.jl`, replace the type loop with behavior variety in one model:
```julia
            behaviorsets = [(), (BrownBerg(),), (Brumley(chemotactic_precision = 0),)]
            for B1 in behaviorsets, B2 in behaviorsets
                space = ContinuousSpace(ntuple(_ -> 100, D); periodic=true)
                model = StandardABM(Microbe{D}, space, 0.1)
                motility = RunTumble([30.0], 0.67, 0.0)
                add_agent!(model; motility, behaviors = B1)
                add_agent!(position(model[1]), model; motility, behaviors = B2)
```
(rest of the loop body unchanged).

- [ ] **Step 7: Run tests**

`julia_restart`; `using Pkg; Pkg.test()`
Expected: PASS including both regression testsets. Agents.jl names `adata` columns after the function name; if the local accessors in `run_scenario` produce different column names than the reference (e.g. mangled closure names), move them to top level and pass `model` through a global `Ref`. If a regression value differs, compare the operation order against the pre-refactor file (`git show 733ef48:src/chemotaxis/<file>.jl`) — any reordering of RNG draws or arithmetic breaks bitwise equality.

- [ ] **Step 8: Commit**

```bash
git add -A src test
git commit -m "feat!: BrownBerg and Brumley as behaviors; add regression harness"
```

---

### Task 6: Celani and Xie as behaviors

**Files:**
- Rewrite: `src/chemotaxis/celani.jl`, `src/chemotaxis/xie.jl`
- Modify: `src/MicrobeAgents.jl` (re-enable their includes)
- Modify: `test/regression.jl`, `test/chemotaxis.jl`, `test/model-creation.jl`

**Interfaces:**
- Consumes: `check_sensing_radius`, `CONV_NOISE` (Task 5), `run_scenario`/`test_against_reference` (Task 5).
- Produces: `Celani(; gain=50.0, memory=1.0, chemotactic_precision=0.0)` with fields `markovian_variables::Vector{Float64}`, `state`; `Xie(; adaptation_time_m=1.29, adaptation_time_z=0.28, gain_forward=2.7, gain_backward=1.6, binding_affinity=0.39, chemotactic_precision=0.0)` with fields `state, state_m, state_z`.

- [ ] **Step 1: Add failing regression + unit tests**

Append inside the regression testset of `test/regression.jl`:
```julia
    test_against_reference("Celani", run_scenario(ref_rt,
        (chemotaxis = Celani(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.26, radius = 0.5), []))
    # accessor names must match the reference column names
    state_m(m) = m.behaviors.chemotaxis.state_m
    state_z(m) = m.behaviors.chemotaxis.state_z
    test_against_reference("Xie", run_scenario(ref_rrf,
        (chemotaxis = Xie(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.26, radius = 0.5), [state_m, state_z]))
```

Append to `test/model-creation.jl` (inside the top-level testset):
```julia
    @testset "Celani initialization at steady state" begin
        for D in 1:3
            C = 2.0
            concentration_field(microbe, model) = C
            chemo = GenericChemoattractant{D}(; concentration_field)
            s = ContinuousSpace(ones(SVector{D}))
            model = StandardABM(Microbe{D}, s, 1.0; properties = Dict(:chemoattractant => chemo))
            add_agent!(model; motility = RunTumble([30.0], 0.67, 0.1), behaviors = (Celani(),))
            c = model[1].behaviors[1]
            λ = 1 / c.memory
            @test c.state == 0.0
            @test c.markovian_variables == [C/λ, C/λ^2, 2C/λ^3]
        end
    end
```
Restore the Celani testset in `test/chemotaxis.jl` with the same mechanical conversion as Task 5 Step 6 (`behaviors = (Celani(gain=5, memory=1),)` etc.), and add to the zero-radius testset:
```julia
        m2 = add_agent!(model; motility, behaviors = (Celani(), Xie()))
        run!(model, 5)
        @test !isnan(bias(m2, model))
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = (Xie(chemotactic_precision = 1),))
```
Add a Xie unit test to `test/chemotaxis.jl`:
```julia
    @testset "Xie forward/backward gains" begin
        space = ContinuousSpace((100.0, 100.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        motility = RunReverseFlick([0.0], Inf, [0.0], Inf)
        add_agent!(model; motility, behaviors = (Xie(),))
        m = model[1]
        x = m.behaviors[1]
        x.state = 0.5
        @test bias(m, model) == 1 + x.gain_forward * 0.5
        m.motility.current_state = 3
        @test bias(m, model) == 1 + x.gain_backward * 0.5
    end
```

- [ ] **Step 2: Run to verify failure**

`julia_restart`; `using MicrobeAgents, Test; include("test/regression.jl")`
Expected: `UndefVarError: Celani`.

- [ ] **Step 3: Rewrite `src/chemotaxis/celani.jl`**

```julia
export Celani

"""
    Celani(; gain=50, memory=1, chemotactic_precision=0)
Chemotaxis behavior using the response kernel from 'Celani and Vergassola (2010) PNAS',
extracted from experiments on E. coli.
Optional sensing noise follows the Berg-Purcell formula, scaled by
`chemotactic_precision` as in 'Brumley et al. (2019) PNAS'; it requires a
microbe `radius > 0`.
Internal markovian variables are initialized at steady state with the local
concentration when the microbe is added.

Parameters:
- `gain = 50`
- `memory = 1` s
- `chemotactic_precision = 0`
Internal state: `state`, `markovian_variables`.
"""
@kwdef mutable struct Celani
    gain::Float64 = 50.0
    memory::Float64 = 1.0
    chemotactic_precision::Float64 = 0.0
    markovian_variables::Vector{Float64} = zeros(3)
    state::Float64 = 0.0
end

function initialize!(b::Celani, microbe, model)
    check_sensing_radius(b, b.chemotactic_precision, microbe)
    W = b.markovian_variables
    λ = 1 / b.memory
    M = concentration(microbe, model)
    W[1] = M / λ
    W[2] = W[1] / λ
    W[3] = 2W[2] / λ
    return nothing
end

function affect!(b::Celani, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = chemoattractant_diffusivity(microbe, model)
    c = concentration(microbe, model)
    a = radius(microbe)
    Π = b.chemotactic_precision
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt)) # noise (Berg-Purcell)
    M = rand(abmrng(model), Normal(c, σ)) # measurement
    λ = 1 / b.memory
    W = b.markovian_variables
    W[1] += (-λ * W[1] + M) * Δt
    W[2] += (-λ * W[2] + W[1]) * Δt
    W[3] += (-λ * W[3] + 2 * W[2]) * Δt
    b.state = λ^2 * (W[2] - λ * W[3] / 2)
    return nothing
end

bias(b::Celani, microbe, model) = 1 - b.gain * b.state
```

- [ ] **Step 4: Rewrite `src/chemotaxis/xie.jl`**

```julia
export Xie

"""
    Xie(; adaptation_time_m=1.29, adaptation_time_z=0.28,
        gain_forward=2.7, gain_backward=1.6, binding_affinity=0.39,
        chemotactic_precision=0)
Chemotaxis behavior adapted from 'Xie et al. (2019) Biophys J', based on the
response function of V. alginolyticus. With a 4-state motility
(`RunReverseFlick`), `gain_backward` applies in the backward run (state 3)
and `gain_forward` otherwise.
Optional Berg-Purcell sensing noise (`chemotactic_precision`) requires a
microbe `radius > 0`.

Parameters:
- `adaptation_time_m = 1.29` s
- `adaptation_time_z = 0.28` s
- `gain_forward = 2.7` 1/s
- `gain_backward = 1.6` 1/s
- `binding_affinity = 0.39` μM
- `chemotactic_precision = 0`
Internal state: `state`, `state_m`, `state_z`.
"""
@kwdef mutable struct Xie
    adaptation_time_m::Float64 = 1.29
    adaptation_time_z::Float64 = 0.28
    gain_forward::Float64 = 2.7
    gain_backward::Float64 = 1.6
    binding_affinity::Float64 = 0.39
    chemotactic_precision::Float64 = 0.0
    state::Float64 = 0.0
    state_m::Float64 = 0.0
    state_z::Float64 = 0.0
end

initialize!(b::Xie, microbe, model) =
    check_sensing_radius(b, b.chemotactic_precision, microbe)

function affect!(b::Xie, microbe::AbstractMicrobe, model)
    Δt = abmtimestep(model)
    Dc = chemoattractant_diffusivity(microbe, model)
    c = concentration(microbe, model)
    K = b.binding_affinity
    a = radius(microbe)
    Π = b.chemotactic_precision
    # noisy concentration measurement with Berg-Purcell formula
    σ = iszero(Π) ? 0.0 : CONV_NOISE * Π * sqrt(3 * c / (5 * π * Dc * a * Δt))
    M = max(rand(abmrng(model), Normal(c, σ)), zero(c))
    ϕ = log(1.0 + M / K)
    τ_m = b.adaptation_time_m
    τ_z = b.adaptation_time_z
    a₀ = (τ_m * τ_z) / (τ_m - τ_z)
    m = b.state_m
    z = b.state_z
    m += (ϕ - m / τ_m) * Δt
    z += (ϕ - z / τ_z) * Δt
    b.state_m = m
    b.state_z = z
    b.state = a₀ * (m / τ_m - z / τ_z)
    return nothing
end

function bias(b::Xie, microbe::AbstractMicrobe{D,N}, model) where {D,N}
    backward = N == 4 && state(motilepattern(microbe)) == 3
    β = backward ? b.gain_backward : b.gain_forward
    return 1 + β * b.state
end
```

Re-enable `include("chemotaxis/celani.jl")` and `include("chemotaxis/xie.jl")` in `src/MicrobeAgents.jl`.

- [ ] **Step 5: Run tests**

`julia_restart`; `using Pkg; Pkg.test()`
Expected: PASS, including Celani and Xie regression.

- [ ] **Step 6: Commit**

```bash
git add -A src test
git commit -m "feat!: Celani and Xie as behaviors"
```

---

### Task 7: SonMenolascina decomposition

**Files:**
- Rewrite: `src/chemotaxis/son-menolascina.jl`
- Modify: `src/MicrobeAgents.jl` (re-enable include)
- Modify: `test/regression.jl`, `test/chemotaxis.jl`

**Interfaces:**
- Consumes: `BrownBerg` (Task 5), hooks.
- Produces: `Chemokinesis(; threshold=0.05, factor=1.3)` (mutable, field `on`); `SpeedDependentTurnRate(; eta=-0.55, ζ=-0.35, θ=1.0, vT=18.88)`; `SpeedDependentFlick(; p0=0.055, amplitude=0.72, steepness=0.25, v_half=36.0)`; `SonMenolascina(; gain, receptor_binding_constant, memory, eta, ζ, θ, vT, threshold, factor)` returning `NamedTuple` `(chemokinesis, chemotaxis, turnrate, flick)` — chemokinesis **first**, because the original model updated chemokinesis before BrownBerg read the velocity.

- [ ] **Step 1: Add failing tests**

Append in `test/regression.jl`'s testset:
```julia
    test_against_reference("SonMenolascina", run_scenario(
        () -> RunReverseFlick([30.0], 0.5, [30.0], 0.5),
        SonMenolascina(),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
```
Restore the SonMenolascina testset in `test/chemotaxis.jl`, converted: `add_agent!(model; motility, behaviors = SonMenolascina(gain=600, memory=1))`, `bias(model[i].behaviors.chemotaxis, model[i], model)` in place of `bias(model[i])` (the old `bias` excluded the speed-dependent factor), and keep `@test speed(model[1]) == 20*1.3`. Add:
```julia
    @testset "Composable SonMenolascina pieces" begin
        bs = SonMenolascina()
        @test keys(bs) == (:chemokinesis, :chemotaxis, :turnrate, :flick)
        extended = (; bs..., extra = (m, model) -> nothing)
        @test length(extended) == 5
        # flick hook only acts when leaving the backward run of a 4-state motility
        space = ContinuousSpace((100.0, 100.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        add_agent!(model; motility = RunReverseFlick([30.0], 1.0, [30.0], 1.0),
            behaviors = (flick = SpeedDependentFlick(),))
        m = model[1]
        w = [0.0, 0.0, 0.0, 1.0]
        m.motility.current_state = 3
        transition_weights!(w, m.behaviors.flick, m, model)
        p = 0.055 + 0.72 / (1 + exp(-0.25 * (30.0 - 36.0)))
        @test w ≈ [0.0, 1 - p, 0.0, p]
        w = [0.0, 0.0, 0.0, 1.0]
        m.motility.current_state = 1
        transition_weights!(w, m.behaviors.flick, m, model)
        @test w == [0.0, 0.0, 0.0, 1.0]
    end
```

- [ ] **Step 2: Run to verify failure**

`julia_restart`; `using MicrobeAgents, Test; include("test/regression.jl")`
Expected: `UndefVarError: SonMenolascina`.

- [ ] **Step 3: Rewrite `src/chemotaxis/son-menolascina.jl`**

```julia
export SonMenolascina, Chemokinesis, SpeedDependentTurnRate, SpeedDependentFlick

"""
    Chemokinesis(; threshold=0.05, factor=1.3)
Multiply the microbe speed by `factor` while the local concentration is
at least `threshold` (μM).
"""
@kwdef mutable struct Chemokinesis
    threshold::Float64 = 0.05
    factor::Float64 = 1.3
    on::Bool = false
end

function affect!(b::Chemokinesis, microbe::AbstractMicrobe, model)
    b.on = concentration(microbe, model) >= b.threshold
    return nothing
end
speed_factor(b::Chemokinesis, microbe) = b.on ? b.factor : 1.0

"""
    SpeedDependentTurnRate(; eta=-0.55, ζ=-0.35, θ=1.0, vT=18.88)
Bias the switching rate by `f(v)/f(0)` with
`f(v) = 1 / (eta / (1 + exp(ζ*(v - vT))) + θ)` and `v = speed(microbe)`
(from 'Son, Menolascina and Stocker (2016) PNAS').
"""
@kwdef struct SpeedDependentTurnRate
    eta::Float64 = -0.55 # s
    ζ::Float64 = -0.35 # s/μm
    θ::Float64 = 1.0 # s
    vT::Float64 = 18.88 # μm/s
end

_turnrate(v, b::SpeedDependentTurnRate) = 1 / (b.eta / (1 + exp(b.ζ * (v - b.vT))) + b.θ)
function bias(b::SpeedDependentTurnRate, microbe, model)
    v = speed(microbe)
    return _turnrate(v, b) / _turnrate(zero(v), b)
end

"""
    SpeedDependentFlick(; p0=0.055, amplitude=0.72, steepness=0.25, v_half=36.0)
With a 4-state motility (`RunReverseFlick`), make the probability of flicking
after a backward run depend on speed:
`p = p0 + amplitude / (1 + exp(-steepness*(speed - v_half)))`;
otherwise the run ends in a reversal.
"""
@kwdef struct SpeedDependentFlick
    p0::Float64 = 0.055
    amplitude::Float64 = 0.72
    steepness::Float64 = 0.25 # s/μm
    v_half::Float64 = 36.0 # μm/s
end

flick_probability(b::SpeedDependentFlick, microbe) =
    b.p0 + b.amplitude / (1 + exp(-b.steepness * (speed(microbe) - b.v_half)))

function transition_weights!(w, b::SpeedDependentFlick, microbe::AbstractMicrobe{D,N}, model) where {D,N}
    if N == 4 && state(motilepattern(microbe)) == 3 # backward run
        p = flick_probability(b, microbe)
        w[2] = 1 - p
        w[4] = p
    end
    return nothing
end

"""
    SonMenolascina(; gain=660, receptor_binding_constant=100, memory=1,
        eta=-0.55, ζ=-0.35, θ=1.0, vT=18.88, threshold=0.05, factor=1.3)
Behaviors of the chemotaxis model from 'Son, Menolascina and Stocker (2016) PNAS':
returns the `NamedTuple`
`(chemokinesis = Chemokinesis(...), chemotaxis = BrownBerg(...),
turnrate = SpeedDependentTurnRate(...), flick = SpeedDependentFlick())`,
to be passed as `behaviors` (use with a `RunReverseFlick` motility).
"""
function SonMenolascina(;
    gain = 660.0, receptor_binding_constant = 100.0, memory = 1.0,
    eta = -0.55, ζ = -0.35, θ = 1.0, vT = 18.88,
    threshold = 0.05, factor = 1.3,
)
    (
        chemokinesis = Chemokinesis(; threshold, factor),
        chemotaxis = BrownBerg(; gain, receptor_binding_constant, memory),
        turnrate = SpeedDependentTurnRate(; eta, ζ, θ, vT),
        flick = SpeedDependentFlick(),
    )
end
```
Re-enable `include("chemotaxis/son-menolascina.jl")`.

- [ ] **Step 4: Run tests**

`julia_restart`; `using Pkg; Pkg.test()`
Expected: PASS. If only the SonMenolascina regression differs at a late time point, check whether the first mismatching row follows a flick decision: the old code updated `ProbabilityWeights.sum` incrementally across calls, so a last-ulp difference in the sum can in principle flip a sample. If (and only if) that is the cause, document it in the test with a comment and compare that scenario up to the step before the divergence; any earlier mismatch is a real bug.

- [ ] **Step 5: Commit**

```bash
git add -A src test
git commit -m "feat!: split SonMenolascina into reusable behaviors"
```

---

### Task 8: Migrate examples and docs

**Files:**
- Modify: `examples/Chemotaxis/2_celani_gauss2D.jl`, `3_xie_response-function.jl`, `4_response_functions.jl`, `5_drift_exponential.jl`, `examples/Pathfinder/1_randomwalk.jl`, `2_chemotaxis.jl`
- Create: `docs/src/behaviors.md`
- Modify: `docs/src/introduction.md`, `docs/src/api.md`, `docs/make.jl`
- (Leave `examples/Encounters/sphere3D.jl` alone: already stale, not built.)

**Interfaces:**
- Consumes: full public API from Tasks 2–7.

- [ ] **Step 1: Migrate examples**

Apply these conversions (read each file fully first; keep narrative text consistent with the code):
- `StandardABM(<Model>{D[,N]}, ...)` → `StandardABM(Microbe{D}, ...)`; `Union{BrownBerg{3},Celani{3}}` → `Microbe{3}`; `add_agent!(BrownBerg{3}, model; ...)` → `add_agent!(model; ...)`.
- Model parameters move into the behavior: e.g. `add_agent!(model; motility, chemotactic_precision=6.0, rotational_diffusivity=0.1)` with Celani → `add_agent!(model; motility, rotational_diffusivity=0.1, radius=0.5, behaviors=(chemotaxis=Celani(chemotactic_precision=6.0),))`. Wherever a model had noise > 0 and relied on the old default `radius = 0.5`, pass `radius = 0.5` explicitly; where it relied on old default `rotational_diffusivity`, pass the old value explicitly (BrownBerg/Brumley/SonMenolascina 0.035; Celani/Xie 0.26).
- `adata = [bias]` → `tumblebias(m) = bias(m, model)` and `adata = [tumblebias]`; column `:bias` → `:tumblebias` in `Analysis.adf_to_matrix` calls.
- Xie example: `adata = [tumblebias, xie_m, xie_z]` with `xie_m(m) = m.behaviors.chemotaxis.state_m`, `xie_z(m) = m.behaviors.chemotaxis.state_z`; update the column names in the `adf_to_matrix` calls; update the prose mentioning `state_m`/`state_z` to say they are fields of the `Xie` behavior.
- `5_drift_exponential.jl`: `adata = [position, velocity, :chemotactic_precision]` → `precision(m) = m.behaviors.chemotaxis.chemotactic_precision` and `adata = [position, velocity, precision]`; update the corresponding column name downstream.
- Pathfinder examples: `BrownBerg{2,2}` → `Microbe{2}`; `1_randomwalk.jl` uses no chemoattractant so add no behaviors; `2_chemotaxis.jl` adds `behaviors = (chemotaxis = BrownBerg(),)` to each `add_agent!` (with explicit `rotational_diffusivity` as already done there).

- [ ] **Step 2: Write `docs/src/behaviors.md`**

```markdown
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

Each agent receives its own copy of the behaviors, so their internal state
is never shared. Plain functions and `Behavior` objects are not copied.

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
bounce!(microbe, model) = position(microbe)[1] > 90 && (microbe.vel = -microbe.vel)
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

## Sharing measurements

`concentration`, `gradient`, `time_derivative` and `chemoattractant_diffusivity`
are evaluated at most once per microbe per step, however many behaviors ask
for them. For other expensive quantities, define a sensor behavior placed
before the behaviors that use it; it stores the value in its own field, and
the others read it by name or through `findbehavior`.

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
```

- [ ] **Step 3: Update `docs/src/introduction.md`, `docs/src/api.md`, `docs/make.jl`**

- `introduction.md` lines ~90–100: replace "All the other subtypes of `AbstractMicrobe` work in a similar way..." and the `@docs BrownBerg ...` block with a short paragraph: "Chemotaxis and other behaviors are attached to microbes through the `behaviors` keyword of `add_agent!`, e.g. `add_agent!(model; motility, behaviors = (BrownBerg(),))`; see [Behaviors](behaviors.md)."
- `introduction.md` stepping section (~lines 118–145): the third bullet becomes "call the `affect!` hook of each behavior (e.g. chemotaxis or other user-defined behavior)"; replace the paragraph about overwriting subroutines with: "Custom behaviors are normally added as behaviors (see [Behaviors](behaviors.md)); the subroutines are exported and can still be reordered or replaced in a custom `agent_step!`."
- `api.md`: add `Microbe` at the top of the Microbes list (keep `state`, which now documents only `state(::Motility)`).
- `docs/make.jl` `pages`: insert `"Behaviors" => "behaviors.md",` after `"Introduction" => "introduction.md",`.
- `docs/src/index.md`: add "Son, Menolascina & Stocker, PNAS 2016" to the list of chemotaxis models.

- [ ] **Step 4: Build docs**

julia-mcp with `env_path` = `/home/riccardo/.julia/dev/MicrobeAgents/.claude/worktrees/behaviors-design/docs`:
```julia
const WT = "/home/riccardo/.julia/dev/MicrobeAgents/.claude/worktrees/behaviors-design"
using Pkg; Pkg.develop(PackageSpec(path = WT)); Pkg.instantiate()
include(joinpath(WT, "docs", "make.jl")) # make.jl cd's into docs/ itself
```
Expected: build completes; no `@docs` missing-docstring errors; all Literate examples execute. Then revert any changes Pkg made to `docs/Project.toml`/`docs/Manifest.toml` if they are tracked (`git status docs`; `git checkout -- docs/Project.toml` if modified) and do not commit `docs/build`.

- [ ] **Step 5: Commit**

```bash
git add -A examples docs/src docs/make.jl
git commit -m "docs: migrate examples and document composable behaviors"
```

---

### Task 9: Final verification

**Files:** none new.

- [ ] **Step 1: Full test suite**

`julia_restart`; `using Pkg; Pkg.test()` — Expected: all pass.

- [ ] **Step 2: Leftover references**

`grep -rnE "chemotaxis!|:affect!|state\(model\[|BrownBerg\{|Celani\{|Brumley\{|Xie\{|SonMenolascina\{" src test examples docs/src` — Expected: no matches outside `test/fixtures/generate_reference.jl`.

- [ ] **Step 3: Allocation check on the hot path**

julia-mcp:
```julia
using MicrobeAgents, Random
space = ContinuousSpace((1000.0, 1000.0))
model = StandardABM(Microbe{2}, space, 0.1; container = Vector)
for _ in 1:100
    add_agent!(model; motility = RunTumble([30.0], 1.0, Isotropic(2)), behaviors = SonMenolascina())
end
run!(model, 10)
@time run!(model, 1000)
```
Expected: allocations dominated by Agents.jl bookkeeping; compare with the same script on `Microbe{2}` without behaviors — the behaviors run should not allocate per agent-step beyond motile switches (≲ a few allocations per switch). Report numbers; investigate if allocations scale with agents × steps.

- [ ] **Step 4: Commit any fixes**

```bash
git add -A
git commit -m "chore: final cleanups for composable behaviors"
```
(skip if nothing changed)
