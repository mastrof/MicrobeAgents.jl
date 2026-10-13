using MicrobeAgents, Test, Random
using Distributions
using LinearAlgebra: norm, ⋅

@testset "Motility" begin
    @test Run([30.0], 0.0) isa MotileState
    @test Tumble(Isotropic2D, 0.0) isa MotileState
    @test Reverse(0.0) isa MotileState
    @test Flick(0.0) isa MotileState
    @test Stop(1.0) isa MotileState
    @test Isotropic(2) === Uniform(-π, π)
    @test Isotropic(3) === Spherical(0, π)

    motile_states = (Run([30.0], 1.0), Tumble(Isotropic2D, 0.1))
    rates = [(1=>2, 1.0), (2=>1, 0.9), (2=>2, 0.1)]
    m = Motility(motile_states, rates...)
    @test m isa Motility{2}
    w = transition_weights(m)
    @test w[1] == [0.0, 1.0]
    @test w[2] == [0.9, 0.1]

    rt = RunTumble(Normal(50, 5.0), 1.0, Isotropic2D; tumble_duration=0.13)
    @test rt isa Motility{2}
    @test transition_weights(rt, 1) == [0.0, 1.0]
    @test transition_weights(rt, 2) == [1.0, 0.0]
    @test duration(states(rt)[1]) == 1.0
    @test duration(states(rt)[2]) == 0.13
    @test speed(states(rt)[1]) == Normal(50, 5.0)
    @test speed(states(rt)[2]) == [0.0]

    rr = RunReverse([40.0], 1.0, [25.0], 0.5)
    @test rr isa Motility{4}
    @test transition_weights(rr, 1) == [0.0, 1.0, 0.0, 0.0]
    @test transition_weights(rr, 2) == [0.0, 0.0, 1.0, 0.0]
    @test transition_weights(rr, 3) == [0.0, 0.0, 0.0, 1.0]
    @test transition_weights(rr, 4) == [1.0, 0.0, 0.0, 0.0]
    @test duration(states(rr)[1]) == 1.0
    @test duration(states(rr)[2]) == 0.0
    @test duration(states(rr)[3]) == 0.5
    @test duration(states(rr)[4]) == 0.0
    @test speed(states(rr)[1]) == [40.0]
    @test speed(states(rr)[2]) == [0.0]
    @test speed(states(rr)[3]) == [25.0]
    @test speed(states(rr)[4]) == [0.0]

    rrf = RunReverseFlick([40.0], 1.0, [25.0], 0.5)
    @test rrf isa Motility{4}
    for i in 1:4
        @test transition_weights(rrf, i) == [j == mod1(i + 1, 4) for j in 1:4]
    end
    @test duration.(states(rrf)) == [1.0, 0.0, 0.5, 0.0]
    rs = RunStop([55.0, 92.0], 1.0, 2.0)
    @test rs isa Motility{2}
    @test transition_weights(rs, 1) == [0.0, 1.0]
    @test transition_weights(rs, 2) == [1.0, 0.0]
    @test duration.(states(rs)) == [1.0, 2.0]
    @test speed(states(rs)[2]) == [0.0]

    # only runs can be biased by behaviors
    @test biased.(states(rrf)) == [true, false, true, false]
    @test biased.(states(rs)) == [true, false]

    @testset "Motile state of a microbe" begin
        model = StandardABM(Microbe{2}, ContinuousSpace((10.0, 10.0)), 0.1; rng = Xoshiro(1))
        m = add_agent!(model; motility = RunReverse([40.0], 1.0, [25.0], 0.5))
        @test motilestate(m) === states(motilepattern(m))[1]
        @test speed(m) == 40.0
        update_motilestate!(motilepattern(m), 3)
        @test state(motilepattern(m)) == 3 && motilestate(m) === states(motilepattern(m))[3]
        update_speed!(m, model)
        @test speed(m) == 25.0
    end

    @testset "turn! samples the reorientation angle" begin
        for D in 2:3, θ in (π/6, π/2, 3π/4)
            model = StandardABM(Microbe{D}, ContinuousSpace(fill(10.0, SVector{D})), 0.1; rng = Xoshiro(2))
            m = add_agent!(model; motility = RunTumble([30.0], 1.0, (θ,)))
            update_motilestate!(motilepattern(m), 2) # tumble
            @test all(1:20) do _
                e = direction(m)
                turn!(m, model)
                norm(direction(m)) ≈ 1 && acos(clamp(direction(m) ⋅ e, -1, 1)) ≈ θ
            end
        end
        # in 1D a turn is always an inversion
        model = StandardABM(Microbe{1}, ContinuousSpace((10.0,)), 0.1)
        m = add_agent!(model; motility = RunTumble([30.0], 1.0, Isotropic(1)))
        e = direction(m)
        turn!(m, model)
        @test direction(m) == -e
    end
end
