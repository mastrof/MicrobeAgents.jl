using MicrobeAgents, Test, Random
using Distributions
using LinearAlgebra: norm

@testset "Motility" begin
    @test AbstractMotilityOneStep <: AbstractMotility
    @test AbstractMotilityTwoStep <: AbstractMotility
    @test ~(AbstractMotilityOneStep <: AbstractMotilityTwoStep)

    @test fieldnames(RunTumble) == (:speed, :polar, :azimuthal)
    rt = RunTumble()
    # type hierarchy
    @test rt isa AbstractMotilityOneStep
    # default values
    @test rt.speed == (30.0,)
    @test rt.polar == Uniform(-π, π)
    @test rt.azimuthal == Arccos()
    # keywords
    rt = RunTumble(; azimuthal=(π,), speed=Normal(5,0.1), polar=[-0.1, 0.1])
    @test rt.speed == Normal(5, 0.1)
    @test rt.polar == [-0.1, 0.1]
    @test rt.azimuthal == (π,)
    # base constructor without keywords
    rt2 = RunTumble(Normal(5,0.1), [-0.1, 0.1], (π,))
    @test rt2.speed == rt.speed && rt2.polar == rt.polar && rt2.azimuthal == rt.azimuthal

    @test fieldnames(RunReverse) == (
        :speed_forward, :polar_forward, :azimuthal_forward,
        :speed_backward, :polar_backward, :azimuthal_backward,
        :motile_state
    )
    # default values
    rr = RunReverse()
    @test rr isa AbstractMotilityTwoStep
    @test rr.speed_forward == (30.0,)
    @test rr.polar_forward == (π,)
    @test rr.azimuthal_forward == Arccos()
    @test rr.speed_backward == rr.speed_forward
    @test rr.polar_backward == rr.polar_forward
    @test rr.azimuthal_backward == rr.azimuthal_forward
    @test rr.motile_state.state == Forward
    # field overload
    @test rr.state == Forward
    # keywords
    rr = RunReverse(; azimuthal_backward = (-π/4,0,π/4))
    @test rr.azimuthal_backward == (-π/4,0,π/4)
    # backward distributions follow forward if unspecified
    rr = RunReverse(; speed_forward = (45,))
    @test rr.speed_backward == rr.speed_forward == (45,)

    @test fieldnames(RunReverseFlick) == fieldnames(RunReverse)
    # default values
    rrf = RunReverseFlick()
    @test rrf isa AbstractMotilityTwoStep
    @test rrf.speed_forward == (30.0,)
    @test rrf.polar_forward == (π,)
    @test rrf.azimuthal_forward == Arccos()
    @test rrf.speed_backward == rrf.speed_forward
    @test rrf.polar_backward == (-π/2, π/2)
    @test rrf.azimuthal_backward == rrf.azimuthal_forward
    @test rrf.motile_state.state == Forward
    # field overload
    @test rrf.state == Forward
    # keywords
    rrf = RunReverseFlick(; azimuthal_backward = (-π/4,0,π/4))
    @test rrf.azimuthal_backward == (-π/4,0,π/4)
    # polar distributions are independent in run reverse flick
    rrf = RunReverseFlick(; speed_forward=(45,), polar_forward=(3π,))
    @test rrf.speed_backward == rrf.speed_forward == (45,)
    @test rrf.polar_forward == (3π,) && rrf.polar_backward == (-π/2,π/2)

    @testset "Motile state" begin
        @test instances(TwoState) == (Forward, Backward)
        # call without argument defaults to Forward
        @test TwoState() == Forward
        # bitwise-not switches between Forward and Backward
        @test ~Forward == Backward
        @test ~Backward == Forward
        # MotileState automatically calls TwoState() if unspecified
        ms₁ = MotileState()
        ms₂ = MotileState(Forward)
        @test ms₁.state == ms₂.state

        rt = RunTumble()
        rr = RunReverse()
        # on a one-step motility nothing happens
        @test_nowarn switch!(rt)
        # on a two-step motility the motile state is inverted
        switch!(rr)
        @test rr.state == Backward
    end

    @testset "Random velocities and angles" begin
        rt = RunTumble()
        @test rand_speed(rt) == 30.0
        rt = RunTumble(speed = Normal(50,5))
        # custom rng
        u₁ = rand(Xoshiro(1234), Normal(50,5))
        u₂ = rand_speed(Xoshiro(1234), rt)
        @test u₁ == u₂

        rr = RunReverse(speed_forward=[45], speed_backward=[35]) # initialized to Forward
        @test rand_speed(rr) == 45
        switch!(rr) # now Backward
        @test rand_speed(rr) == 35
        # custom rng
        rr = RunReverse(speed_forward = Normal(10,1))
        u₁ = rand(Random.seed!(1), Normal(10,1))
        u₂ = rand_speed(Random.seed!(1), rr)
        @test u₁ == u₂

        rt = RunTumble()
        for D in 1:3
            v = rand_vel(D, rt)
            @test length(v) == D && norm(v) ≈ 30.0
        end
        # custom rng
        rt = RunTumble(speed = Uniform(30,40))
        for D in 1:3
            v = rand_vel(Random.seed!(1), D, rt)
            Random.seed!(1)
            # generate and discard D values corresponding to the velocity projection
            rand(D)
            # the D+1-th value is extracted from the speed distribution
            u = rand(Uniform(30,40))
            @test norm(v) ≈ u
        end

        rt = RunTumble(speed=[25], polar=[0.1], azimuthal=[0.2])
        U, θ, ϕ = rand(rt)
        @test U==25 && θ==0.1 && ϕ==0.2
        rt = RunTumble()
        rng = Random.seed!(Xoshiro(1), 27)
        U, θ, ϕ = rand(rng, rt)
        rng = Random.seed!(Xoshiro(1), 27)
        U′, θ′, ϕ′ = rand(rng,(30.0,)), rand(rng,Uniform(-π,π)), rand(rng,Arccos())
        @test U==U′ && θ==θ′ && ϕ==ϕ′

        rr = RunReverseFlick() # initialized to Forward
        rng = Random.seed!(Xoshiro(1), 42)
        U, θ, ϕ = rand(rng, rr)
        rng = Random.seed!(Xoshiro(1), 42)
        U′, θ′, ϕ′ = rand(rng,(30.0,)), rand(rng,(π,)), rand(rng,Arccos())
        @test U==U′ && θ==θ′ && ϕ==ϕ′
        switch!(rr) # now Backward
        rng = Random.seed!(Xoshiro(1), 23)
        U, θ, ϕ = rand(rng, rr)
        rng = Random.seed!(Xoshiro(1), 23)
        U′, θ′, ϕ′ = rand(rng,(30.0,)), rand(rng,(-π/2,π/2)), rand(rng,Arccos())
        @test U==U′ && θ==θ′ && ϕ==ϕ′
    end
end