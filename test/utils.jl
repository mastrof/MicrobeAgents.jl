using MicrobeAgents, Test
using StaticArrays: SVector
using LinearAlgebra
using Random

@testset "Utility functions" begin
    for D in 1:3
        extent = ones(SVector{D})
        space = ContinuousSpace(extent)
        model = StandardABM(Microbe{D}, space)
        v = random_velocity(model)
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
end
