using MicrobeAgents, Test, Random
using Distributions
using LinearAlgebra: norm

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

    rrf = RunReverseFlick([40.0], 1.0, [25.0], 1.0)
    @test rr isa Motility{4}
    rs = RunStop([55.0, 92.0], 1.0, 2.0)
    @test rs isa Motility{2}
end
