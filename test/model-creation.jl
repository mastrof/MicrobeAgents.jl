using MicrobeAgents, Test
using LinearAlgebra: norm
using Random

@testset "Model creation" begin
    for D in 1:3
        timestep = 1
        space = ContinuousSpace(ones(SVector{D}))
        model = StandardABM(Microbe{D}, space, timestep)
        @test model isa StandardABM
        @test Set(keys(abmproperties(model))) == Set((
            :timestep,
            :chemoattractant,
            :field_cache
        ))
    end

    @testset "Base Microbe type" begin
        for D in 1:3
            timestep = 1
            space = ContinuousSpace(ones(SVector{D}))
            model = StandardABM(Microbe{D}, space, timestep; rng=Xoshiro(123))
            # add agent with default constructor
            # random pos and vel, random speed from motility pattern
            motility = RunTumble([25.0, 35.0], 1.3, 0.1)
            add_agent!(model; motility)
            rng = Xoshiro(123)
            pos = rand(rng, SVector{D})
            vel = random_velocity(rng, D)
            #speed = random_speed(rng, RunTumble())
            spd = rand(rng, speed(motility))
            @test model[1] isa Microbe{D}
            @test model[1] isa Microbe{D,2}
            @test position(model[1]) == pos
            @test direction(model[1]) == vel
            @test speed(model[1]) == spd
            @test radius(model[1]) == 0.0
            @test rotational_diffusivity(model[1]) == 0.0
            # add agent with predefined position
            pos = SVector{D}(i/2D for i in 1:D)
            add_agent!(pos, model; motility)
            @test model[2].pos == pos
        end
    end
end
