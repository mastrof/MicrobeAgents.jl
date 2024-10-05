using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

@testset "Model stepping" begin
    @testset "Instantaneous turns" begin
        for D in 1:3, container in (Dict, Vector)
            rng = Xoshiro(68)
            dt = 1
            extent = fill(300.0, SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt; rng, container)
            pos = extent ./ 2
            m1 = RunTumble(Inf, [30.0], Isotropic(D)) # infinite run
            vel1 = random_velocity(model)
            speed1 = rand(rng, speed(m1))
            add_agent!(pos, model; vel=vel1, speed=speed1, motility=m1)
            m2 = RunReverse(0, [10.0], 0, [10.0]) # switching every step
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

            # customize microbe affect! function
            # decreases microbe state value by D at each step
            affect!(microbe::Microbe{D}, model) where D = (microbe.state -= D)
            properties = Dict(:affect! => affect!)
            model = StandardABM(Microbe{D}, space, dt; container, properties)
            motility = RunTumble(1.0, [30.0], Isotropic(D), 0.0)
            add_agent!(model; motility)
            run!(model, 1)
            @test model[1].state == -D

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
            m1 = RunTumble(Inf, [30.0], Isotropic(D), 0.1) # infinite run
            vel1 = random_velocity(model)
            speed1 = rand(rng, speed(m1))
            add_agent!(pos, model; vel=vel1, speed=speed1, motility=m1)
            m2 = RunReverse(0, [10.0], 0, [10.0], 0.1, 0.1) # switching every step
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

            # customize microbe affect! function
            # decreases microbe state value by D at each step
            affect!(microbe::Microbe{D}, model) where D = (microbe.state -= D)
            properties = Dict(:affect! => affect!)
            model = StandardABM(Microbe{D}, space, dt; container, properties)
            motility = RunTumble(1.0, [30.0], Isotropic(D), 0.0)
            add_agent!(model; motility)
            run!(model, 1)
            @test model[1].state == -D

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
