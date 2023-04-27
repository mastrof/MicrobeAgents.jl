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

    @testset "BrownBerg" begin
        for D in 1:3
            m = BrownBerg{D}()
            @test typeof(m) == BrownBerg{D}
            @test m isa AbstractMicrobe{D}
            @test fieldnames(BrownBerg{D}) == (
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :rotational_diffusivity,
                :radius, :state, :gain,
                :receptor_binding_constant, :memory
            )
            @test m.pos isa NTuple{D,Float64}
            @test m.vel isa NTuple{D,Float64}
            @test m.speed == 30.0
            @test m.motility isa RunTumble

            model = StandardABM(BrownBerg{D}, ntuple(_->100,D), 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m, model) == m.turn_rate*exp(-m.gain*m.state)
            @test m.state == 0
        end
    end
    @testset "Brumley" begin
        for D in 1:3
            m = Brumley{D}()
            @test typeof(m) == Brumley{D}
            @test m isa AbstractMicrobe{D}
            @test fieldnames(Brumley{D}) == (
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :rotational_diffusivity,
                :radius, :state, :memory,
                :gain_receptor, :gain,
                :chemotactic_precision
            )
            @test m.pos isa NTuple{D,Float64}
            @test m.vel isa NTuple{D,Float64}
            @test m.speed == 46.5
            @test m.motility isa RunReverseFlick

            model = StandardABM(Brumley{D}, ntuple(_->100,D), 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m, model) == (1+exp(-m.gain*m.state))*m.turn_rate/2
            @test m.state == 0.0
        end
    end
    @testset "Celani" begin
        for D in 1:3
            m = Celani{D}()
            @test typeof(m) == Celani{D}
            @test m isa AbstractMicrobe{D}
            @test fieldnames(Celani{D}) == (
                :id, :pos, :motility, :vel, :speed,
                :turn_rate, :markovian_variables, :state,
                :rotational_diffusivity, :gain, :memory,
                :radius, :chemotactic_precision
            )
            @test m.pos isa NTuple{D,Float64}
            @test m.vel isa NTuple{D,Float64}
            @test m.speed == 30.0
            @test m.motility isa RunTumble
            
            model = StandardABM(Celani{D}, ntuple(_->100,D), 0.1)
            add_agent!(model)
            m = model[1]
            @test turnrate(m,model) == m.turn_rate*(1-m.gain*m.state)
            affect!(m, model)
            # concentration_field is null, so state should be unchanged
            @test m.markovian_variables == zeros(3)
            @test m.state == 1.0
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