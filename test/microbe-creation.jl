using MicrobeAgents, Test
using LinearAlgebra: norm

@testset "Microbe creation" begin
    for D in 1:3
        m = Microbe{D}()
        # type hierarchy
        @test typeof(m) == Microbe{D}
        @test m isa AbstractMicrobe{D}
        # correct fields for all Ds
        @test fieldnames(Microbe{D}) == (
            :id, :pos, :motility, :vel, :speed, :turn_rate,
            :rotational_diffusivity, :radius, :state
        )
        # pos and vel sizes match D
        @test m.pos isa NTuple{D,Float64}
        @test m.vel isa NTuple{D,Float64}
        @test m.speed == 30.0
        # default values
        @test m.motility isa RunTumble
        @test m.turn_rate == 1.0
        for key in fieldnames(Microbe{D})[end:-1:end-2]
            @test getfield(m, key) == 0.0
        end
    end
    # keyword tests
    motility = RunReverse(speed_forward=[42.0], speed_backward=[60.1])
    m = Microbe{3}(1, (2,9,0), radius=0.5; motility) # id and pos don't require keyword
    @test m.id == 1 && m.pos == (2.0,9.0,0.0) && m.radius == 0.5
    @test m.motility isa RunReverse
    @test m.motility.speed_forward == [42.0]
    @test m.motility.speed_backward == [60.1]
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

    @testset "Celani" begin
        for D in 1:3
            m = Celani{D}()
            @test typeof(m) == Celani{D}
            @test m isa AbstractMicrobe{D}
            @test fieldnames(Celani{D}) == (
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :state, :rotational_diffusivity,
                :gain, :memory, :radius, :chemotactic_precision
            )
            @test m.pos isa NTuple{D,Float64}
            @test m.vel isa NTuple{D,Float64}
            @test m.speed == 30.0
            @test m.motility isa RunTumble
            # can use arbitrary 2nd argument since model is not called
            @test turnrate(m,0) == m.turn_rate*m.state[4]
        end
    end
    @testset "Xie" begin
        for D in 1:3
            m = Xie{D}()
            @test typeof(m) == Xie{D}
            @test m isa AbstractMicrobe{D}
            @test fieldnames(Xie{D}) == (
                :id, :pos, :motility, :vel, :speed,
                :turn_rate_forward, :turn_rate_backward,
                :rotational_diffusivity, :radius,
                :state, :state_m, :state_z,
                :adaptation_time_m, :adaptation_time_z,
                :gain_forward, :gain_backward,
                :binding_affinity, :chemotactic_precision
            )
            @test m.pos isa NTuple{D,Float64}
            @test m.vel isa NTuple{D,Float64}
            @test m.speed == 46.5
            @test m.motility isa RunReverseFlick
        end
    end
end