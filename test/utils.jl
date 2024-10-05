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

    @testset "Distances" begin
        for D in 1:3
            MicrobeTypes = [Microbe{D}, Celani{D}, BrownBerg{D}, Brumley{D}]
            for T1 in MicrobeTypes, T2 in MicrobeTypes
                space = ContinuousSpace(ntuple(_ -> 100, D); periodic=true)
                model = StandardABM(Union{T1,T2}, space, 0.1)
                motility = RunTumble([30.0], 0.67, 0.0)
                add_agent!(T1, model; motility)
                add_agent!(position(model[1]), T2, model; motility)
                delta = SVector{D}(randn(D) .* 5)
                walk!(model[2], delta, model)
                @test distance(model[1], model[2], model) ≈ norm(delta)
                v = distancevector(model[1], model[2], model)
                @test norm(v) ≈ norm(delta)
            end
        end
    end
end
