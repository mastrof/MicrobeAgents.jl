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
            :affect!
        ))
    end

    @testset "Base Microbe type" begin
        for D in 1:3
            timestep = 1
            space = ContinuousSpace(ones(SVector{D}))
            model = StandardABM(Microbe{D}, space, timestep; rng=Xoshiro(123))
            # add agent with default constructor
            # random pos and vel, random speed from motility pattern
            motility = RunTumble(1.3, [25.0, 35.0], 0.1)
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
            @test state(model[1]) == 0.0
            # add agent with predefined position
            pos = SVector{D}(i/2D for i in 1:D)
            add_agent!(pos, model; motility)
            @test model[2].pos == pos
        end
    end

    @testset "Chemotactic Microbe types" begin
        for T in [BrownBerg, Celani, Brumley], D in 1:3
            @testset "$(T{D})" begin
                timestep = 1
                space = ContinuousSpace(ones(SVector{D}))
                model = StandardABM(T{D}, space, timestep; rng=Xoshiro(987))
                rng = Xoshiro(987)
                motility = RunTumble(1.0, [30.0], 0.0)
                add_agent!(model; motility)
                m = model[1]
                @test m isa T{D}
                @test position(m) == rand(rng, SVector{D})
                @test direction(m) == random_velocity(rng, D)
                @test speed(m) == rand(rng, speed(motilepattern(m)))
                @test issubset(
                    (:id, :pos, :vel, :speed, :motility,
                    :rotational_diffusivity, :radius, :state),
                    fieldnames(T)
                )

                if T == Celani
                    # when no concentration field is set, markovian variables are zero
                    @test m.markovian_variables == zeros(3)

                    # initialize a new model with non-zero concentration field
                    C = 2.0
                    concentration_field(pos, model) = C
                    chemo = GenericChemoattractant{D,Float64}(;concentration_field)
                    properties = Dict(:chemoattractant => chemo)
                    s = ContinuousSpace(ones(SVector{D}))
                    model_c = StandardABM(Celani{D}, s, 1.0; properties)
                    motility = RunTumble(0.67, [30.0], 0.1)
                    add_agent!(model_c; motility)
                    m = model_c[1]
                    λ = 1 / m.memory
                    @test m.state == 0.0
                    @test m.markovian_variables == [C/λ, C/λ^2, 2C/λ^3]
                end
            end
        end
    end
end
