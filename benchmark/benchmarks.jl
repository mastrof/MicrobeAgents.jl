using BenchmarkTools
using MicrobeAgents
using Random

# Benchmark suite for MicrobeAgents.jl, following the PkgBenchmark/AirspeedVelocity
# convention: this file must define a `SUITE::BenchmarkGroup`.
#
# Run and compare revisions locally with AirspeedVelocity.jl:
#   benchpkg --rev=dirty,main            # working tree vs main
#   benchpkg --rev=v0.6.0,v0.7.0         # across releases
#   benchpkg table --rev=v0.6.0,v0.7.0   # (after benchpkg) print results
# Pull requests are compared automatically by .github/workflows/benchmark.yml.
#
# Only the package itself and BenchmarkTools are used, so that the suite can
# run against old revisions with minimal friction.

const SUITE = BenchmarkGroup()

const DIMENSIONS = (1, 2, 3)
const DT = 0.1
const BOXSIZE = 100.0
const SEED = 1234

# microbe types and the motility pattern (number of states) they require
const NSTATES = Dict(
    Microbe => 2, BrownBerg => 2, Celani => 2, Xie => 4,
    Brumley => 4, SonMenolascina => 4,
)

function make_motility(T, D)
    if NSTATES[T] == 2
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

function make_model(T, D, N; chemotaxis = false)
    space = ContinuousSpace(ntuple(_ -> BOXSIZE, D); periodic = true)
    properties = chemotaxis ? Dict(:chemoattractant => linear_chemoattractant(D)) : Dict()
    model = StandardABM(
        T{D,NSTATES[T]}, space, DT;
        container = Vector, rng = Xoshiro(SEED), properties,
    )
    for _ in 1:N
        add_agent!(model; motility = make_motility(T, D))
    end
    model
end

# `evals = 1` plus a fresh model in `setup` keeps samples independent.
# `run!(model, nsteps)` is amortised over several steps to beat timer noise.
function stepping_benchmark(T, D, N; chemotaxis = false, nsteps = 100)
    @benchmarkable(
        run!(model, $nsteps),
        setup = (model = make_model($T, $D, $N; chemotaxis = $chemotaxis)),
        evals = 1,
    )
end

# Single microbe, no chemoattractant: measures bare stepping overhead
SUITE["stepping"] = BenchmarkGroup()
for T in keys(NSTATES)
    group = SUITE["stepping"]["$(nameof(T))"] = BenchmarkGroup()
    for D in DIMENSIONS
        group["$(D)D"] = stepping_benchmark(T, D, 1)
    end
end

# Chemotactic microbes in a non-trivial field
SUITE["chemotaxis"] = BenchmarkGroup()
for T in (BrownBerg, Brumley, Celani, Xie, SonMenolascina)
    group = SUITE["chemotaxis"]["$(nameof(T))"] = BenchmarkGroup()
    for D in DIMENSIONS
        group["$(D)D"] = stepping_benchmark(T, D, 1; chemotaxis = true)
    end
end

# Population scaling in 3D
SUITE["population"] = BenchmarkGroup()
for T in (Microbe, BrownBerg)
    group = SUITE["population"]["$(nameof(T))"] = BenchmarkGroup()
    for N in (10, 1000)
        group["N=$N"] =
            stepping_benchmark(T, 3, N; chemotaxis = T != Microbe, nsteps = 10)
    end
end

# Model construction
SUITE["setup"] = BenchmarkGroup()
for D in DIMENSIONS
    SUITE["setup"]["add_agent! $(D)D"] = @benchmarkable(
        add_agent!(model; motility = motility),
        setup = (
            model = make_model(Microbe, $D, 0);
            motility = make_motility(Microbe, $D)
        ),
        evals = 1,
    )
end
