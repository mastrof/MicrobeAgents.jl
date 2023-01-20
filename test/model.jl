using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

@testset "Agent Based Model" begin
    @testset "Model creation" begin
        for D in 1:3
            timestep = 1
            extent = ntuple(_ -> 300, D)
            model = ABM(Microbe{D}, extent, timestep)
            @test Set(keys(model.properties)) == Set((
                :t, :timestep,
                :concentration_field,
                :concentration_gradient,
                :concentration_time_derivative,
                :compound_diffusivity,
                :update!
            ))

            # add agents with default constructor, random pos and vel
            rng = MersenneTwister(1234)
            model = ABM(Microbe{D}, extent, timestep; rng)
            add_agent!(model)
            rng = MersenneTwister(1234)
            pos = Tuple(rand(rng,D) .* extent)
            vel = rand_vel(rng, D)
            speed = rand_speed(rng, RunTumble())
            @test Set(keys(model.agents)) == Set((1,))
            @test model[1] isa Microbe{D}
            @test model[1].pos == pos
            @test model[1].vel == vel
            @test model[1].speed == speed
            @test model[1].motility isa RunTumble
            m₀ = Microbe{D}()
            for key in setdiff(fieldnames(Microbe{D}), (:id,:pos,:vel,:motility))
                @test getfield(model[1], key) == getfield(m₀, key)
            end

            # add agents with keyword arguments
            model = ABM(Microbe{D}, extent, timestep)
            motility = RunReverse(speed_backward=[24.0], motile_state=MotileState(Backward))
            add_agent!(model; turn_rate=0.55, motility)
            @test model[1].turn_rate == 0.55
            @test model[1].motility isa RunReverse
            @test norm(model[1].vel) ≈ 1.0
            @test model[1].speed ≈ 24.0
            vel = rand_vel(D)
            speed = rand_speed(RunTumble())
            add_agent!(model; vel)
            @test model[2].vel == vel
            @test model[2].speed == speed

            # add agent to given position
            model = ABM(Microbe{D}, extent, timestep)
            pos = extent ./ 2
            add_agent!(pos, model)
            @test model[1] isa Microbe{D}
            @test model[1].pos == pos
            pos = extent ./ 3
            microbe = Microbe{D}(25) # id=25
            add_agent!(microbe, pos, model)
            @test model[25].pos == pos
            # conserve agent's own position
            pos = extent ./ 6
            microbe = Microbe{D}(2, pos)
            add_agent_pos!(microbe, model)
            @test model[2].pos == pos
        end
    end

    @testset "Stepping" begin
        for D in 1:3
            dt = 1
            extent = ntuple(_ -> 300, D)
            model = ABM(Microbe{D}, extent, dt)
            pos = extent ./ 2
            vel1 = rand_vel(D)
            speed1 = rand_speed(RunTumble())
            add_agent!(pos, model; vel=vel1, speed=speed1, turn_rate=0)
            vel2 = rand_vel(D)
            speed2 = rand_speed(RunReverse())
            add_agent!(pos, model; vel=vel2, speed=speed2, turn_rate=Inf, motility=RunReverse())
            run!(model) # performs 1 microbe_step! and 1 model.update! step
            # x₁ = x₀ + vΔt
            @test all(model[1].pos .≈ pos .+ vel1.*speed1.*dt)
            @test all(model[2].pos .≈ pos .+ vel2.*speed2.*dt)
            # v is the same for the agent with zero turn rate
            @test all(model[1].vel .≈ vel1)
            # v is changed for the other agent
            @test ~all(model[2].vel .≈ vel2)
            # and since it turned its motile state changed from Forward to Backward
            @test model[2].motility.state == Backward
            # model time counter was increased to 1
            @test model.t == 1
            # perform other 5 steps
            run!(model, 5)
            @test model.t == 6

            # customize microbe affect! function
            # decreases microbe state value by D at each step
            MicrobeAgents.affect!(microbe::Microbe{D}, model) = (microbe.state-=D)
            model = ABM(Microbe{D}, extent, dt)
            add_agent!(model)
            run!(model)
            @test model[1].state == -D
            # restore default affect! function
            MicrobeAgents.affect!(microbe::Microbe{D}, model) = MicrobeAgents._affect!(microbe, model)
            run!(model)
            # now the state has not changed from previous step
            @test model[1].state == -D

            # customize model.update! function
            model = ABM(Microbe{D}, extent, dt)
            tick_more!(model::ABM) = (model.t += 3)
            model → tick_more! # now model.t increases by 4 (+1 +3) at each step
            run!(model, 6)
            @test model.t == 6*4
            decrease!(model::ABM) = (model.t -= 1)
            model = ABM(Microbe{D}, extent, dt)
            # chain arbitrary number of functions
            model → tick_more! → decrease! → decrease! → decrease!
            # now we should be back at only +1 per step
            run!(model, 20)
            @test model.t == 20
        end
    end
end