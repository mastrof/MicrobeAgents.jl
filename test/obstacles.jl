using MicrobeAgents, Test
using LinearAlgebra: norm, dot
using Random

const MA = MicrobeAgents

# finite-difference normal (generic fallback), to check analytic normals
fd_normal(p, x::SVector{D}) where {D} =
    invoke(MA.surface_normal, Tuple{MA.AbstractObstacle{D},SVector{D,Float64}}, p, x)

@testset "Obstacles" begin
    @testset "Primitives" begin
        s = ObstacleSphere((0.0, 0.0), 2.0)
        @test s isa ObstacleSphere{2}
        @test sdf(s, SVector(3.0, 4.0)) ≈ 3.0
        @test sdf(s, SVector(1.0, 0.0)) ≈ -1.0
        @test MA.surface_normal(s, SVector(3.0, 4.0)) ≈ SVector(0.6, 0.8)
        @test MA.boundingradius(s) == 2.0
        @test sdf(ObstacleSphere([1.0], 0.5), SVector(2.0)) ≈ 0.5 # 1D
        @test sdf(ObstacleSphere((0, 0, 0), 1), SVector(0.0, 0.0, 3.0)) ≈ 2.0 # 3D

        h = ObstacleHalfSpace((0.0, 0.0), (0.0, 2.0))
        @test sdf(h, SVector(5.0, 3.0)) ≈ 3.0
        @test MA.surface_normal(h, SVector(5.0, 3.0)) ≈ SVector(0.0, 1.0)
        @test MA.boundingradius(h) == Inf

        b = ObstacleBox((0.0, 0.0), (2.0, 4.0))
        @test sdf(b, SVector(3.0, 2.0)) ≈ 1.0
        @test MA.surface_normal(b, SVector(3.0, 2.0)) ≈ SVector(1.0, 0.0)
        @test sdf(b, SVector(1.0, 2.0)) ≈ -1.0
        @test sdf(b, SVector(3.0, 6.0)) ≈ sqrt(5)
        @test MA.surface_normal(b, SVector(3.0, 6.0)) ≈ SVector(1.0, 2.0) / sqrt(5)
        @test MA.boundingradius(b) ≈ sqrt(5)
        @test_throws ArgumentError ObstacleBox((2.0, 0.0), (0.0, 1.0))

        c = ObstacleCapsule((0.0, 0.0), (4.0, 0.0), 1.0)
        @test sdf(c, SVector(2.0, 3.0)) ≈ 2.0
        @test MA.surface_normal(c, SVector(2.0, 3.0)) ≈ SVector(0.0, 1.0)
        @test sdf(c, SVector(6.0, 0.0)) ≈ 1.0
        @test MA.boundingradius(c) ≈ 3.0
        @test MA.center(c) ≈ SVector(2.0, 0.0)

        f = ObstacleSDF(x -> norm(x) - 1, (0.0, 0.0), 1.0)
        @test sdf(f, SVector(2.0, 0.0)) ≈ 1.0
        @test MA.surface_normal(f, SVector(2.0, 0.0)) ≈ SVector(1.0, 0.0) atol = 1e-6
        @test !MA.isconvex(f)
        @test MA.isconvex(ObstacleSDF(x -> norm(x) - 1, (0.0, 0.0), 1.0; convex = true))
        @test all(MA.isconvex, (s, h, b, c))

        # analytic normals agree with the finite-difference fallback
        rng = Xoshiro(3)
        for D in 1:3
            ps = (ObstacleSphere(fill(5.0, D), 2.0), ObstacleBox(fill(3.0, D), fill(6.0, D)),
                ObstacleCapsule(fill(3.0, D), fill(7.0, D), 1.0), ObstacleHalfSpace(fill(5.0, D), ones(D)))
            for p in ps, _ in 1:200
                x = SVector{D}(10 .* rand(rng, D))
                abs(sdf(p, x)) < 1e-3 && continue
                sdf(p, x) < 0 && continue # gradient ill-defined inside (medial axes, 1D capsule)
                @test MA.surface_normal(p, x) ≈ fd_normal(p, x) atol = 1e-4
            end
        end
    end
end
