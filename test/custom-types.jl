using MicrobeAgents, Test, Random
using LinearAlgebra: norm
using MicrobeAgents.Agents: @agent, ContinuousAgent

# user-defined microbe type carrying behaviors through a `behaviors(m)` method
@agent struct Swimmer{D,N}(ContinuousAgent{D,Float64}) <: AbstractMicrobe{D,N}
    speed::Float64
    motility::Motility{N}
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
    energy::Float64 = 10.0
end
struct Fatigue end
MicrobeAgents.affect!(::Fatigue, m::Swimmer, model) = (m.energy -= 1; nothing)
MicrobeAgents.speed_factor(::Fatigue, m::Swimmer) = m.energy / 10
MicrobeAgents.behaviors(::Swimmer) = (Fatigue(),)

@testset "Custom microbe types" begin
    for D in 1:3
        space = ContinuousSpace(fill(1000.0, SVector{D}))
        model = StandardABM(Swimmer{D}, space, 0.1; rng = Xoshiro(1))
        motility = RunTumble([30.0], Inf, Isotropic(D))
        m = add_agent!(fill(500.0, SVector{D}), model; motility)
        @test m isa Swimmer{D,2}
        @test m.energy == 10.0 && speed(m) == 30.0
        x0 = position(m)
        run!(model, 4)
        @test m.energy == 6.0
        @test speed(m) ≈ 18.0
        @test findbehavior(m, Fatigue) === Fatigue()
        # displacement uses the speed at each step (30, 27, 24, 21) × dt
        @test norm(position(m) - x0) ≈ (30 + 27 + 24 + 21) * 0.1
    end
end
