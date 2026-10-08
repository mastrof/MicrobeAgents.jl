using BenchmarkTools
using MicrobeAgents
using Random

# Benchmark suite for MicrobeAgents.jl, following the PkgBenchmark/AirspeedVelocity
# convention: this file must define a `SUITE::BenchmarkGroup`.
#
# Run and compare revisions locally with AirspeedVelocity.jl:
#   benchpkg --rev=dirty,main            # working tree vs main
#   benchpkg --rev=v1.0.0,dirty          # across revisions (v1+ only)
#   benchpkg table --rev=v1.0.0,dirty    # (after benchpkg) print results
# Pull requests are compared automatically by .github/workflows/benchmark.yml.
#
# Only the package itself and BenchmarkTools are used.
# The suite targets the behaviors API (MicrobeAgents v1+): chemotactic models
# are `behaviors` of the single `Microbe` type, so it cannot run on older revisions.

const SUITE = BenchmarkGroup()

const DIMENSIONS = (1, 2, 3)
const DT = 0.1
const BOXSIZE = 100.0
const SEED = 1234
const RADIUS = 0.5 # finite radius, required by behaviors with sensing noise

# model name => (number of motile states, constructor of fresh behaviors)
const MODELS = Dict(
    "Microbe" => (2, () -> ()),
    "BrownBerg" => (2, () -> (chemotaxis = BrownBerg(),)),
    "Celani" => (2, () -> (chemotaxis = Celani(),)),
    "Xie" => (4, () -> (chemotaxis = Xie(),)),
    "Brumley" => (4, () -> (chemotaxis = Brumley(),)),
    "SonMenolascina" => (4, () -> SonMenolascina()),
)
nstates(name) = MODELS[name][1]
make_behaviors(name) = MODELS[name][2]()

function make_motility(name, D)
    if nstates(name) == 2
        RunTumble(
            run_speed = [30.0], run_duration = 1.0, angle = Isotropic(D)
        )
    else
        RunReverseFlick(
            run_speed_forward = [30.0], run_duration_forward = 1.0,
            run_speed_backward = [30.0], run_duration_backward = 1.0,
        )
    end
end

# constant-gradient field along x, so that chemotactic models do real work
function linear_chemoattractant(D)
    c(m, model) = 1.0 + 0.01 * m.pos[1]
    g(m, model) = SVector{D,Float64}(ntuple(i -> i == 1 ? 0.01 : 0.0, D))
    r(m, model) = 0.01 * m.vel[1] * m.speed
    GenericChemoattractant{D}(
        concentration_field = c,
        concentration_gradient = g,
        concentration_ramp = r,
    )
end

function make_model(name, D, N; chemotaxis = false)
    space = ContinuousSpace(ntuple(_ -> BOXSIZE, D); periodic = true)
    properties = chemotaxis ? Dict(:chemoattractant => linear_chemoattractant(D)) : Dict()
    # concrete agent type, as recommended in the Behaviors docs
    T = Microbe{D,nstates(name),typeof(make_behaviors(name))}
    model = StandardABM(
        T, space, DT;
        container = Vector, rng = Xoshiro(SEED), properties,
    )
    for _ in 1:N
        add_agent!(
            model; motility = make_motility(name, D),
            behaviors = make_behaviors(name), radius = RADIUS,
        )
    end
    model
end

# `evals = 1` plus a fresh model in `setup` keeps samples independent.
# `run!(model, nsteps)` is amortised over several steps to beat timer noise.
function stepping_benchmark(name, D, N; chemotaxis = false, nsteps = 100)
    @benchmarkable(
        run!(model, $nsteps),
        setup = (model = make_model($name, $D, $N; chemotaxis = $chemotaxis)),
        evals = 1,
    )
end

# Single microbe, no chemoattractant: measures bare stepping overhead
SUITE["stepping"] = BenchmarkGroup()
for name in keys(MODELS)
    group = SUITE["stepping"][name] = BenchmarkGroup()
    for D in DIMENSIONS
        group["$(D)D"] = stepping_benchmark(name, D, 1)
    end
end

# Chemotactic microbes in a non-trivial field
SUITE["chemotaxis"] = BenchmarkGroup()
for name in ("BrownBerg", "Brumley", "Celani", "Xie", "SonMenolascina")
    group = SUITE["chemotaxis"][name] = BenchmarkGroup()
    for D in DIMENSIONS
        group["$(D)D"] = stepping_benchmark(name, D, 1; chemotaxis = true)
    end
end

# Population scaling in 3D
SUITE["population"] = BenchmarkGroup()
for name in ("Microbe", "BrownBerg")
    group = SUITE["population"][name] = BenchmarkGroup()
    for N in (10, 1000)
        group["N=$N"] =
            stepping_benchmark(name, 3, N; chemotaxis = name != "Microbe", nsteps = 10)
    end
end

# Model construction; `add_agent!` also runs `initialize!` on the behaviors
SUITE["setup"] = BenchmarkGroup()
for name in ("Microbe", "Celani", "SonMenolascina"), D in DIMENSIONS
    SUITE["setup"]["add_agent! $name $(D)D"] = @benchmarkable(
        add_agent!(
            model; motility = motility, behaviors = behaviors, radius = RADIUS,
        ),
        setup = (
            model = make_model($name, $D, 0);
            motility = make_motility($name, $D);
            behaviors = make_behaviors($name)
        ),
        evals = 1,
    )
end
