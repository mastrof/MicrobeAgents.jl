using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

# custom per-step logic: type and method must be defined at top level
mutable struct Countdown
    s::Float64
end
MicrobeAgents.affect!(c::Countdown, microbe::Microbe{D}, model) where {D} = (c.s -= D)

@testset "Model stepping" begin
    @testset "Translation and turns (turn duration = $τ_turn)" for τ_turn in (0.0, 0.1)
        for D in 1:3
            rng = Xoshiro(68)
            dt = 1
            extent = fill(300.0, SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt; rng)
            pos = extent ./ 2
            m1 = RunTumble(;
                run_duration=Inf, run_speed=[30.0],
                angle=Isotropic(D), tumble_duration=τ_turn
            ) # infinite run
            vel1 = random_velocity(model)
            speed1 = rand(rng, speed(m1))
            add_agent!(pos, model; vel=vel1, speed=speed1, motility=m1)
            m2 = RunReverse(;
                run_duration_forward=0, run_speed_forward=[10.0],
                run_duration_backward=0, run_speed_backward=[10.0],
                reverse_duration=τ_turn,
            ) # switching every step
            vel2 = random_velocity(model)
            speed2 = rand(rng, speed(m2))
            add_agent!(pos, model; vel=vel2, speed=speed2, motility=m2)
            run!(model, 1)
            # x₁ = x₀ + vΔt
            @test position(model[1]) ≈ @. pos + vel1 * speed1 * dt
            @test position(model[2]) ≈ @. pos + vel2 * speed2 * dt
            # v is the same for the agent with zero turn rate
            @test velocity(model[1]) ≈ vel1 .* speed1
            # instantaneous reversals complete within the step; finite ones leave the microbe still
            @test velocity(model[2]) ≈ (τ_turn == 0 ? -vel2 .* speed2 : zero(vel2))
        end
    end

    @testset "Stop" begin
        space = ContinuousSpace((100.0, 100.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        motility = RunStop([30.0], Inf, Inf)
        m = add_agent!(SVector(50.0, 50.0), model; motility)
        update_motilestate!(motilepattern(m), 2)
        update_speed!(m, model)
        run!(model, 10)
        @test speed(m) == 0
        @test position(m) == SVector(50.0, 50.0)
    end

    @testset "Custom per-step logic as a behavior" begin
        for D in 1:3
            space = ContinuousSpace(fill(300.0, SVector{D}))
            model = StandardABM(Microbe{D}, space, 1)
            motility = RunTumble(run_duration=1.0, run_speed=[30.0], angle=Isotropic(D))
            add_agent!(model; motility, behaviors = (Countdown(0.0),))
            run!(model, 2)
            @test model[1].behaviors[1].s == -2D
        end
    end
end
