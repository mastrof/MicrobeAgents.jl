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
                model = StandardABM(S, $T, dt; container=Vector),
                S = ContinuousSpace(SVector{$D}(10.0 for _ in $D)),
                dt = 0.1
            )
        )
    end
end
==#

version = Pkg.project().version
filename = "benchmarks_$(version).md"
open(filename, "w") do io
    @printf io "# MicrobeAgents.jl benchmarks - %s \n\n" version
end

## Setup
microbe_types = [Microbe, BrownBerg, Brumley, Celani]
dimensions = [1, 2, 3]

##
open(filename, "a") do io
    @printf io "## Simple stepping\n"
    @printf io "%-16s | %-16s | %-16s | %-16s\n" "Type" "Time (μs)" "Allocs" "Memory"
    @printf io " --- | --- | --- | --- \n"
    for T in microbe_types
        for D in dimensions
            space = ContinuousSpace(ntuple(_ -> 10.0, D))
            dt = 0.1
            if T == Brumley
                model = StandardABM(T{D,4}, space, dt; container=Vector)
                motility = RunReverseFlick(0.0, [10.0], 0.0, [10.0])
            else
                model = StandardABM(T{D,2}, space, dt; container=Vector)
                motility = RunTumble(0.0, [10.0], 0)
            end
            add_agent!(model; motility)
            nsteps = 1000
            b = @benchmark run!($model, $(nsteps))
            type = string(T{D})
            mintime = minimum(b.times)*1e-3 / nsteps
            allocs = b.allocs / nsteps
            memory = b.memory / nsteps
            @printf io "%-16s | %-16.2e | %-16.2e | %-16.2e\n" type mintime allocs memory
        end
    end
end
