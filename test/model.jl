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
            vel = rand_vel(rng, D, RunTumble())
            @test Set(keys(model.agents)) == Set((1,))
            @test model[1] isa Microbe{D}
            @test model[1].pos == pos
            @test model[1].vel == vel
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
            @test norm(model[1].vel) ≈ 24.0

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
        
    end
end