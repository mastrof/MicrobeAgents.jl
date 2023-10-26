using BenchmarkTools
using MicrobeAgents
using Printf
import Pkg

#==
##
const SUITE = BenchmarkGroup()
SUITE["Simple stepping"] = BenchmarkGroup()

##
microbe_types = [Microbe, BrownBerg, Brumley, Celani, Xie]
dimensions = [1, 2, 3]

for D in dimensions
    SUITE["Simple stepping"]["$D dimensions"] = BenchmarkGroup()
    for T in microbe_types
        SUITE["Simple stepping"]["$D dimensions"]["$T"] = @benchmarkable(
            run!(model, 1),
            evals = 10,
            samples = 1000,
            setup = (
                model = UnremovableABM(S, $T, dt),
                S = ContinuousSpace(SVector{$D}(10.0 for _ in $D)),
                dt = 0.1
            )
        )
    end
end
==#

open("benchmarks.md", "w") do io
    @printf io "# MicrobeAgents.jl benchmarks - %s \n\n" Pkg.project().version
end

## Setup
microbe_types = [Microbe, BrownBerg, Brumley, Celani, Xie]
dimensions = [1, 2, 3]

##
open("benchmarks.md", "a") do io
    @printf io "## Simple stepping\n"
    @printf io "%-16s %-16s %-16s %-16s\n" "Type" "Time (μs)" "Allocs" "Memory"
    for T in microbe_types
        for D in dimensions
            space = ContinuousSpace(ntuple(_ -> 10.0, D))
            dt = 0.1
            model = UnremovableABM(T{D}, space, dt)
            b = @benchmark run!($model, $(1))
            @printf io "%-16s %-16.2e %-16d %-16.2e\n" string(T{D}) minimum(b.times)*1e-3 b.allocs b.memory
        end
    end
end
