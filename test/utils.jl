using MicrobeAgents, Test
using LinearAlgebra
using Random

@testset "Utility functions" begin
    for D in 1:3
        v = rand_vel(D)
        @test norm(v) ≈ 1 && length(v) == D
        rng = MersenneTwister()
        v = rand_vel(rng, D)
        @test norm(v) ≈ 1 && length(v) == D
    end

    g₁!(a) = (a .+= 1)
    g₂!(a) = (a .+= 2)
    f₁! = ChainedFunction(g₁!, g₂!)
    f₂! = g₁! → g₂!
    @test f₁! === f₂!
    a₀ = [0,1]; a₁ = [0,1]
    g₁!(a₀); g₂!(a₀)
    f₁!(a₁)
    @test a₀ == a₁ == [3,4]

    g₃!(a) = (a .*= 2)
    f₁! = ChainedFunction(g₁!, ChainedFunction(g₂!, ChainedFunction(g₃!, g₁!)))
    f₂! = g₁! → g₂! → g₃! → g₁!
    @test f₁! === f₂!
    a = [0,1]
    f₁!(a)
    @test a == [7,9]

    @testset "Distances" begin
        L = 20.0
        dt = 1.0
        for D in 1:3
            extent = ntuple(_ -> L, D)

            # periodic
            model = StandardABM(Microbe{D}, extent, dt)
            x₁ = ntuple(i -> i==1 ? 1.0 : 0.0, D)
            x₂ = ntuple(i -> i==1 ? L-1 : 0.0, D)
            add_agent!(x₁, model)
            add_agent!(x₂, model)
            p₁ = ntuple(i -> i==1 ? 2.0 : 0.0, D)
            p₂ = ntuple(i -> i==1 ? L-2 : 0.0, D)
            # microbe-microbe
            @test distance(model[1], model[2], model) ≈ 2
            # microbe-point
            @test distance(model[1], p₁, model) ≈ 1
            @test distance(model[1], p₂, model) ≈ 3
            # point-point
            @test distance(p₁, p₂, model) ≈ 4

            # closed box
            # periodic
            model = StandardABM(Microbe{D}, extent, dt; periodic=false)
            x₁ = ntuple(i -> i==1 ? 1.0 : 0.0, D)
            x₂ = ntuple(i -> i==1 ? L-1 : 0.0, D)
            add_agent!(x₁, model)
            add_agent!(x₂, model)
            p₁ = ntuple(i -> i==1 ? 2.0 : 0.0, D)
            p₂ = ntuple(i -> i==1 ? L-2 : 0.0, D)
            # microbe-microbe
            @test distance(model[1], model[2], model) ≈ L-2
            # microbe-point
            @test distance(model[1], p₁, model) ≈ 1
            @test distance(model[1], p₂, model) ≈ L-3
            # point-point
            @test distance(p₁, p₂, model) ≈ L-4
        end
    end
end