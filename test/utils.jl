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
end
