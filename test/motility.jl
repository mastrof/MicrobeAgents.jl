using MicrobeAgents, Test, Random
using Distributions
using LinearAlgebra: norm

@testset "Motility" begin
    @test MotilityOneStep <: AbstractMotility
    @test MotilityTwoStep <: AbstractMotility
    @test !(MotilityOneStep <: MotilityTwoStep)

    @test fieldnames(MotilityOneStep) == (:speed, :polar, :azimuthal)
    rt = RunTumble()
    # type hierarchy
    @test rt isa MotilityOneStep
    # default values
    @test rt.speed == (30.0,)
    @test rt.polar == Uniform(-π, π)
    @test rt.azimuthal == Arccos()
    # keywords
    rt = RunTumble(; azimuthal=(π,), speed=Normal(5,0.1), polar=[-0.1, 0.1])
    @test rt.speed == Normal(5, 0.1)
    @test rt.polar == [-0.1, 0.1]
    @test rt.azimuthal == (π,)

    @test fieldnames(MotilityTwoStep) == (
        :speed, :polar, :azimuthal,
        :speed_backward, :polar_backward, :azimuthal_backward,
        :motile_state
    )
    # default values
    rr = RunReverse()
    @test rr isa MotilityTwoStep
    @test rr.speed == (30.0,)
    @test rr.polar == (π,)
    @test rr.azimuthal == Arccos()
    @test rr.speed_backward == rr.speed
    @test rr.polar_backward == rr.polar
    @test rr.azimuthal_backward == rr.azimuthal
    @test rr.motile_state.state == Forward
    # field overload
    @test rr.state == Forward
    # keywords
    rr = RunReverse(; azimuthal_backward = (-π/4,0,π/4))
    @test rr.azimuthal_backward == (-π/4,0,π/4)
    # backward distributions follow forward if unspecified
    rr = RunReverse(; speed = (45,))
    @test rr.speed_backward == rr.speed == (45,)

    # default values
    rrf = RunReverseFlick()
    @test rrf isa MotilityTwoStep
    @test rrf.speed == (30.0,)
    @test rrf.polar == (π,)
    @test rrf.azimuthal == Arccos()
    @test rrf.speed_backward == rrf.speed
    @test rrf.polar_backward == (-π/2, π/2)
    @test rrf.azimuthal_backward == rrf.azimuthal
    @test rrf.motile_state.state == Forward
    # field overload
    @test rrf.state == Forward
    # keywords
    rrf = RunReverseFlick(; azimuthal_backward = (-π/4,0,π/4))
    @test rrf.azimuthal_backward == (-π/4,0,π/4)
    # polar distributions are independent in run reverse flick
    rrf = RunReverseFlick(; speed=(45,), polar=(3π,))
    @test rrf.speed_backward == rrf.speed == (45,)
    @test rrf.polar == (3π,) && rrf.polar_backward == (-π/2,π/2)

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
        model = StandardABM(Microbe{2}, ContinuousSpace((1,1)))
        rt = RunTumble()
        add_agent!(model; motility = rt)
        @test random_speed(model[1], model) == 30.0
        # custom rng
        rt = RunTumble(speed = Normal(50,5))
        u₁ = rand(Xoshiro(1234), Normal(50,5))
        model = StandardABM(Microbe{2}, ContinuousSpace((1,1)); rng=Xoshiro(1234))
        # set pos, vel, speed to not trigger the rng
        add_agent!((0,0), model; vel=(0,0), speed=50, motility=rt)
        u₂ = random_speed(model[1], model)
        @test u₁ == u₂

        model = StandardABM(Microbe{2}, ContinuousSpace((1,1)))
        rr = RunReverse(speed=[45], speed_backward=[35]) # initialized to Forward
        add_agent!(model; motility = rr)
        @test random_speed(model[1], model) == 45
        switch!(model[1].motility) # now Backward
        @test random_speed(model[1], model) == 35

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
