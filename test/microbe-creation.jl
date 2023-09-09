using MicrobeAgents, Test
using LinearAlgebra: norm

@testset "Microbe creation" begin
    for D in 1:3
        m = Microbe{D}(; id=1, pos=zero(SVector{D}), vel=zero(SVector{D}), speed=0.0)
        # type hierarchy
        @test typeof(m) == Microbe{D}
        @test m isa AbstractMicrobe{D}
        # correct fields for all Ds
        @test Set(fieldnames(Microbe{D})) == Set((
            :id, :pos, :motility, :vel, :speed, :turn_rate,
            :rotational_diffusivity, :radius, :state
        ))
        # pos and vel sizes match D
        @test m.pos isa SVector{D,Float64}
        @test m.vel isa SVector{D,Float64}
        @test m.speed == 0.0
        # default values
        @test m.motility isa RunTumble
        @test m.turn_rate == 1.0
        for key in fieldnames(Microbe{D})[end:-1:end-2]
            @test getfield(m, key) == 0.0
        end
    end
    # keyword tests
    motility = RunReverse(speed_forward=[42.0], speed_backward=[60.1])
    m = Microbe{3}(; id=1, pos=(2,9,0), vel=(0,0,1), speed=30, radius=0.5, motility)
    @test m.id == 1 && m.pos == SVector(2.0,9.0,0.0) && m.radius == 0.5
    @test m.motility isa RunReverse
    @test m.motility.speed_forward == [42.0]
    @test m.motility.speed_backward == [60.1]
    #==
    # the initial speed should be sampled from the appropriate speed distribution
    # by default it is Forward
    @test norm(m.vel) ≈ 1.0
    @test m.speed ≈ 42.0
    # now initialize from Backward
    motility = RunReverse(speed_forward=[42.0], speed_backward=[60.1],
        motile_state=MotileState(Backward))
    m = Microbe{3}(;motility)
    @test norm(m.vel) ≈ 1
    @test m.speed ≈ 60.1
    ==#

    @testset "BrownBerg" begin
        for D in 1:3
            m = BrownBerg{D}(;id=1,pos=zero(SVector{D}),vel=zero(SVector{D}),speed=0)
            @test typeof(m) == BrownBerg{D}
            @test m isa AbstractMicrobe{D}
            @test Set(fieldnames(BrownBerg{D})) == Set((
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :rotational_diffusivity,
                :radius, :state, :gain,
                :receptor_binding_constant, :memory
            ))
            @test m.pos isa SVector{D,Float64}
            @test m.vel isa SVector{D,Float64}
            @test m.speed == 0.0
            @test m.motility isa RunTumble

            space = ContinuousSpace(ntuple(_->100,D))
            model = StandardABM(BrownBerg{D}, space, 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m, model) == m.turn_rate*exp(-m.gain*m.state)
            @test m.state == 0
        end
    end
    @testset "Brumley" begin
        for D in 1:3
            m = Brumley{D}(;id=1,pos=zero(SVector{D}),vel=zero(SVector{D}),speed=0)
            @test typeof(m) == Brumley{D}
            @test m isa AbstractMicrobe{D}
            @test Set(fieldnames(Brumley{D})) == Set((
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :rotational_diffusivity,
                :radius, :state, :memory,
                :gain_receptor, :gain,
                :chemotactic_precision
            ))
            @test m.pos isa SVector{D,Float64}
            @test m.vel isa SVector{D,Float64}
            @test m.speed == 0
            @test m.motility isa RunReverseFlick

            space = ContinuousSpace(ntuple(_->100,D))
            model = StandardABM(Brumley{D}, space, 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m, model) == (1+exp(-m.gain*m.state))*m.turn_rate/2
            @test m.state == 0.0
        end
    end
    @testset "Celani" begin
        for D in 1:3
            m = Celani{D}(;id=1,pos=zero(SVector{D}),vel=zero(SVector{D}),speed=0)
            @test typeof(m) == Celani{D}
            @test m isa AbstractMicrobe{D}
            @test Set(fieldnames(Celani{D})) == Set((
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :markovian_variables, :state,
                :rotational_diffusivity, :gain, :memory,
                :radius, :chemotactic_precision
            ))
            @test m.pos isa SVector{D,Float64}
            @test m.vel isa SVector{D,Float64}
            @test m.speed == 0.0
            @test m.motility isa RunTumble

            space = ContinuousSpace(ntuple(_->100,D))
            model = StandardABM(Celani{D}, space, 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m,model) == m.turn_rate*(1-m.gain*m.state)
            affect!(m, model)
            # concentration_field is null, so state should be unchanged
            @test m.markovian_variables == zeros(3)
            @test m.state == 0.0

            # if the model has a non-zero concentration field
            # markovian variables and internal state should be
            # initialized from the appropriate steady state values
            C = 2.0
            concentration_field(pos, model) = C
            properties = Dict(:concentration_field => concentration_field)
            space = ContinuousSpace(ntuple(_->100,D))
            model = StandardABM(Celani{D}, space, 0.1; properties)
            add_agent!(model)
            m = model[1]
            p = m.pos
            λ = 1/m.memory
            @test m.state == 0.0
            @test m.markovian_variables == [C/λ, C/λ^2, 2C/λ^3]
        end
    end
    @testset "Xie" begin
        for D in 1:3
            m = Xie{D}(;id=1,pos=zero(SVector{D}),vel=zero(SVector{D}),speed=0)
            @test typeof(m) == Xie{D}
            @test m isa AbstractMicrobe{D}
            @test Set(fieldnames(Xie{D})) == Set((
                :id, :pos, :motility, :vel, :speed,
                :turn_rate_forward, :turn_rate_backward,
                :rotational_diffusivity, :radius,
                :state, :state_m, :state_z,
                :adaptation_time_m, :adaptation_time_z,
                :gain_forward, :gain_backward,
                :binding_affinity, :chemotactic_precision
            ))
            @test m.pos isa SVector{D,Float64}
            @test m.vel isa SVector{D,Float64}
            @test m.speed == 0
            @test m.motility isa RunReverseFlick
        end
    end
end
