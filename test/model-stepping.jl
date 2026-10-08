using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

# custom per-step logic: type and method must be defined at top level
mutable struct Countdown
    s::Float64
end
MicrobeAgents.affect!(c::Countdown, microbe::Microbe{D}, model) where {D} = (c.s -= D)

@testset "Model stepping" begin
    @testset "Instantaneous turns" begin
        for D in 1:3, container in (Dict, Vector)
            rng = Xoshiro(68)
            dt = 1
            extent = fill(300.0, SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt; rng, container)
            pos = extent ./ 2
            m1 = RunTumble(
                run_duration=Inf, run_speed=[30.0], angle=Isotropic(D)
            ) # infinite run
            vel1 = random_velocity(model)
            speed1 = rand(rng, speed(m1))
            add_agent!(pos, model; vel=vel1, speed=speed1, motility=m1)
            m2 = RunReverse(
                run_duration_forward=0, run_speed_forward=[10.0],
                run_duration_backward=0, run_speed_backward=[10.0]
            ) # switching every step
            vel2 = random_velocity(model)
            speed2 = rand(rng, speed(m2))
            add_agent!(pos, model; vel=vel2, speed=speed2, motility=m2)
            run!(model, 1) # performs 1 microbe_step!
            # x₁ = x₀ + vΔt
            @test position(model[1]) ≈ @. pos + vel1 * speed1 * dt
            @test position(model[2]) ≈ @. pos + vel2 * speed2 * dt
            # v is the same for the agent with zero turn rate
            @test velocity(model[1]) ≈ vel1 .* speed1
            # v is changed for the other agent
            @test velocity(model[2]) ≈ -vel2 .* speed2

            # custom per-step logic as a behavior
            model = StandardABM(Microbe{D}, space, dt; container)
            motility = RunTumble(
                run_duration=1.0, run_speed=[30.0],
                angle=Isotropic(D), tumble_duration=0.0
            )
            add_agent!(model; motility, behaviors = (Countdown(0.0),))
            run!(model, 1)
            @test model[1].behaviors[1].s == -D

            # customize model_step! function
            properties = Dict(:square_t => [0])
            model_step!(model) = (abmproperties(model)[:square_t][1] = abmtime(model)^2)
            model = StandardABM(Microbe{D}, space, dt; model_step!, container, properties)
            n = 6
            run!(model, n)
            @test abmproperties(model)[:square_t][1] == (n-1)^2
        end
    end
    @testset "Finite turn times" begin
        for D in 1:3, container in (Dict, Vector)
            rng = Xoshiro(68)
            dt = 1
            extent = fill(300.0, SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt; rng, container)
            pos = extent ./ 2
            m1 = RunTumble(;
                run_duration=Inf, run_speed=[30.0],
                angle=Isotropic(D), tumble_duration=0.1
            ) # infinite run
            vel1 = random_velocity(model)
            speed1 = rand(rng, speed(m1))
            add_agent!(pos, model; vel=vel1, speed=speed1, motility=m1)
            m2 = RunReverse(;
                run_duration_forward=0, run_speed_forward=[10.0],
                run_duration_backward=0, run_speed_backward=[10.0],
                reverse_duration = 0.1,

            ) # switching every step
            vel2 = random_velocity(model)
            speed2 = rand(rng, speed(m2))
            add_agent!(pos, model; vel=vel2, speed=speed2, motility=m2)
            run!(model, 1) # performs 1 microbe_step!
            # x₁ = x₀ + vΔt
            @test position(model[1]) ≈ @. pos + vel1 * speed1 * dt
            @test position(model[2]) ≈ @. pos + vel2 * speed2 * dt
            # v is the same for the agent with zero turn rate
            @test velocity(model[1]) ≈ vel1 .* speed1
            # v is changed for the other agent
            @test velocity(model[2]) ≈ zero(vel2)

            # custom per-step logic as a behavior
            model = StandardABM(Microbe{D}, space, dt; container)
            motility = RunTumble(
                run_duration=1.0, run_speed=[30.0],
                angle=Isotropic(D), tumble_duration=0.0
            )
            add_agent!(model; motility, behaviors = (Countdown(0.0),))
            run!(model, 1)
            @test model[1].behaviors[1].s == -D

            # customize model_step! function
            properties = Dict(:square_t => [0])
            model_step!(model) = (abmproperties(model)[:square_t][1] = abmtime(model)^2)
            model = StandardABM(Microbe{D}, space, dt; model_step!, container, properties)
            n = 6
            run!(model, n)
            @test abmproperties(model)[:square_t][1] == (n-1)^2
        end
    end

end
