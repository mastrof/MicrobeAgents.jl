using MicrobeAgents, Test
using LinearAlgebra: norm, dot
using Random
using MicrobeAgents.Agents: allagents

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
        # SVector arguments must not hit the implicit (center, halfwidth) constructor
        @test sdf(ObstacleBox(SVector(1.0, 1.0), SVector(3.0, 5.0)), SVector(4.0, 3.0)) ≈ 1.0
        @test ObstacleHalfSpace(SVector(0.0, 0.0), SVector(0.0, 2.0)).normal ≈ SVector(0.0, 1.0)

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

    @testset "Contact rules" begin
        space = ContinuousSpace((50.0, 50.0))
        model = StandardABM(Microbe{2}, space, 0.1; rng = Xoshiro(1))
        motility = RunTumble(run_speed = [10.0], run_duration = 1.0, angle = Isotropic(2))
        m = add_agent!(SVector(25.0, 25.0), model; motility, vel = SVector(-0.6, 0.8))
        n = SVector(1.0, 0.0)
        s = SVector(-3.0, 4.0)
        c = SurfaceContact{2}(1, 1, SVector(24.5, 25.0), n, s)

        @test contact!(Halt(), m, c, model) == SVector(0.0, 0.0)
        @test direction(m) == SVector(-0.6, 0.8)

        @test contact!(Slide(), m, c, model) ≈ SVector(0.0, 4.0)
        @test direction(m) == SVector(-0.6, 0.8)
        away = SurfaceContact{2}(1, 1, SVector(24.5, 25.0), n, SVector(3.0, 4.0))
        @test contact!(Slide(), m, away, model) ≈ SVector(3.0, 4.0) # not into the surface
        @test contact!(Slide(align = true), m, c, model) ≈ SVector(0.0, 4.0)
        @test direction(m) ≈ SVector(0.0, 1.0)

        m.vel = SVector(-0.6, 0.8)
        @test contact!(Reflect(), m, c, model) ≈ SVector(3.0, 4.0)
        @test direction(m) ≈ SVector(0.6, 0.8)

        for _ in 1:100
            m.vel = SVector(-0.6, 0.8)
            u = contact!(Randomize(), m, c, model)
            @test dot(direction(m), n) ≥ 0
            @test norm(direction(m)) ≈ 1
            @test norm(u) ≈ norm(s)
            @test u ≈ norm(s) * direction(m)
        end

        @test contact!((mic, con, mod) -> 2 .* con.displacement, m, c, model) == 2 .* s
    end

    # brute-force signed distance over all primitives and periodic images
    function brute_sdf(prims, x::SVector{D}, extent, periodic) where {D}
        best = Inf
        images = Iterators.product(ntuple(d -> periodic[d] ? (-extent[d], 0.0, extent[d]) : (0.0,), D)...)
        for p in prims, s in images
            best = min(best, sdf(p, x - SVector{D}(s)))
        end
        best
    end

    @testset "Obstacle sets and nearest" begin
        rng = Xoshiro(1)
        for D in 1:3, periodic in (false, true), (cellsize, margin) in ((nothing, nothing), (2.0, 0.5), (7.0, 0.0))
            L = 100.0
            space = ContinuousSpace(ntuple(_ -> L, D); periodic)
            prims = [ObstacleSphere(rand(rng, SVector{D}) * L, 2 + 3rand(rng)) for _ in 1:60]
            set = Obstacles(space, prims; cellsize, margin)
            @test set.index isa MA.CellIndex
            for _ in 1:500
                x = rand(rng, SVector{D}) * L
                cutoff = rand(rng, (0.3, 3.0, 12.0, 40.0, Inf))
                d, n, i = MA.nearest(set, x, cutoff)
                b = brute_sdf(prims, x, set.extent, set.periodic)
                if b ≤ cutoff
                    @test d ≈ b # includes points inside solids (b < 0)
                    @test 1 ≤ i ≤ length(prims)
                else
                    @test d == Inf && i == 0
                end
            end
        end

        # obstacle straddling a periodic boundary is seen from the other side
        space = ContinuousSpace((50.0, 50.0); periodic = true)
        set = Obstacles(space, [ObstacleSphere((49.0, 25.0), 3.0)])
        d, n, i = MA.nearest(set, SVector(5.0, 25.0), 10.0)
        @test d ≈ 3.0 && n ≈ SVector(1.0, 0.0) && i == 1

        # walls
        space = ContinuousSpace((50.0, 40.0); periodic = (true, false))
        walls = Walls(space)
        @test walls.index isa MA.NoIndex
        @test length(walls.primitives) == 2
        d, n, i = MA.nearest(walls, SVector(10.0, 3.0), Inf)
        @test d ≈ 3.0 && n ≈ SVector(0.0, 1.0) && i == 1
        d, n, i = MA.nearest(walls, SVector(10.0, 39.0), Inf)
        @test d ≈ 1.0 && n ≈ SVector(0.0, -1.0) && i == 2
        @test_throws ArgumentError Walls(ContinuousSpace((50.0, 50.0); periodic = true))
        @test Walls(space; contact = Slide()).contact == Slide()

        # tuple of sets
        sets = (Obstacles(space, [ObstacleSphere((25.0, 20.0), 5.0)]), walls)
        d, n, k, i = MA._nearest(sets, SVector(25.0, 2.0), Inf)
        @test d ≈ 2.0 && k == 2 && i == 1
        d, n, k, i = MA._nearest(sets, SVector(25.0, 12.0), Inf)
        @test d ≈ 3.0 && k == 1 && i == 1

        # invalid inputs
        @test_throws ArgumentError Obstacles(space, ObstacleSphere{2}[])
        mixed = AbstractObstacle{2}[ObstacleSphere((1.0, 1.0), 1.0), ObstacleBox((3.0, 3.0), (4.0, 4.0))]
        @test_throws ArgumentError Obstacles(space, mixed)
        @test_throws ArgumentError Obstacles(space, [ObstacleSphere((1.0, 1.0), 1.0)]; margin = -1.0)
    end

    @testset "Model integration" begin
        space = ContinuousSpace((50.0, 50.0); periodic = false)
        walls = Walls(space)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:obstacles => walls))
        @test abmproperties(model).obstacles === (walls,) # single set wrapped
        sets = (Obstacles(space, [ObstacleSphere((25.0, 25.0), 5.0)]), walls)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:obstacles => sets))
        @test abmproperties(model).obstacles === sets
        other = ContinuousSpace((60.0, 50.0); periodic = false)
        @test_throws ArgumentError StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => (Walls(other),)))
        @test_throws ArgumentError StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => (Walls(ContinuousSpace((50.0, 50.0); periodic = (true, false))),)))
        @test_throws ArgumentError StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => [walls]))

        # nearest_surface: walls are ordered (x lower, x upper, y lower, y upper)
        motility = RunTumble(run_speed = [10.0], run_duration = 1.0, angle = Isotropic(2))
        m = add_agent!(SVector(8.0, 25.0), model; motility, radius = 0.5)
        q = nearest_surface(m, model)
        @test q.distance ≈ 7.5 && q.normal ≈ SVector(1.0, 0.0) && q.set == 2 && q.index == 1
        m2 = add_agent!(SVector(25.0, 33.0), model; motility, radius = 0.5)
        q = nearest_surface(m2, model)
        @test q.distance ≈ 2.5 && q.normal ≈ SVector(0.0, 1.0) && q.set == 1 && q.index == 1
        @test nearest_surface(m, model; cutoff = 5.0).distance == Inf
        empty_model = StandardABM(Microbe{2}, space, 0.1)
        m3 = add_agent!(SVector(10.0, 25.0), empty_model; motility)
        @test nearest_surface(m3, empty_model) == (distance = Inf, normal = SVector(0.0, 0.0), set = 0, index = 0)

        # random_free_position
        model = StandardABM(Microbe{2}, space, 0.1; rng = Xoshiro(2), properties = Dict(:obstacles => sets))
        for _ in 1:200
            x = random_free_position(model; radius = 1.0)
            @test sdf(sets[1].primitives[1], x) ≥ 1.0
            @test all(w -> sdf(w, x) ≥ 1.0, walls.primitives)
        end
        full = (Obstacles(space, [ObstacleBox((-1.0, -1.0), (51.0, 51.0))]),)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:obstacles => full))
        @test_throws ErrorException random_free_position(model; maxtries = 10)
    end

    # one microbe moving straight (no rotational diffusion, never tumbles)
    function straight_model(space, obstacles, pos, vel; speed = 10.0, radius = 0.5, dt = 1.0)
        model = StandardABM(Microbe{length(pos)}, space, dt; rng = Xoshiro(1),
            properties = Dict(:obstacles => obstacles))
        motility = RunTumble(run_speed = [speed], run_duration = Inf, angle = Isotropic(length(pos)))
        m = add_agent!(SVector{length(pos)}(pos), model; motility, vel = SVector{length(pos)}(vel), radius)
        return model, m
    end

    @testset "Contacts in stepping" begin
        space = ContinuousSpace((50.0, 50.0); periodic = false)
        v = SVector(-1.0, 1.0) / sqrt(2)
        model, m = straight_model(space, (Walls(space),), (2.0, 10.0), v)
        run!(model, 1)
        @test m.pos ≈ SVector(0.5, 11.5) # halted at the wall after moving 1.5 along x
        @test direction(m) ≈ v
        run!(model, 1)
        @test m.pos ≈ SVector(0.5, 11.5) # stays while pointing into the wall

        model, m = straight_model(space, (Walls(space; contact = Slide()),), (2.0, 10.0), v)
        run!(model, 1)
        @test m.pos ≈ SVector(0.5, 10.0 + 10 / sqrt(2))
        @test direction(m) ≈ v

        model, m = straight_model(space, (Walls(space; contact = Slide(align = true)),), (2.0, 10.0), v)
        run!(model, 1)
        @test m.pos ≈ SVector(0.5, 10.0 + 10 / sqrt(2))
        @test direction(m) ≈ SVector(0.0, 1.0)

        model, m = straight_model(space, (Walls(space; contact = Reflect()),), (2.0, 10.0), v)
        run!(model, 1)
        @test m.pos ≈ SVector(0.5 + (10 / sqrt(2) - 1.5), 10.0 + 10 / sqrt(2))
        @test direction(m) ≈ SVector(1.0, 1.0) / sqrt(2)

        model, m = straight_model(space, (Walls(space; contact = Randomize()),), (2.0, 10.0), v)
        run!(model, 1)
        @test direction(m)[1] ≥ 0
        @test m.pos[1] ≥ 0.5

        # obstacle across the periodic boundary, approached from the other side
        pspace = ContinuousSpace((50.0, 50.0); periodic = true)
        model, m = straight_model(pspace, (Obstacles(pspace, [ObstacleSphere((49.0, 25.0), 3.0)]),),
            (5.0, 25.0), (-1.0, 0.0))
        run!(model, 1)
        @test m.pos ≈ SVector(2.5, 25.0)

        # 1D walls
        space1 = ContinuousSpace((20.0,); periodic = false)
        model, m = straight_model(space1, (Walls(space1),), (5.0,), (-1.0,))
        run!(model, 1)
        @test m.pos ≈ SVector(0.5)

        # a (moving) microbe placed overlapping an obstacle is pushed out on its
        # first step (overlap counts as a contact: `Halt` stops it there)
        model, m = straight_model(space, (Obstacles(space, [ObstacleSphere((25.0, 25.0), 5.0)]),),
            (27.0, 25.0), (0.0, 1.0))
        run!(model, 1)
        @test m.pos ≈ SVector(30.5, 25.0)
        run!(model, 1) # then leaves along the tangent
        @test m.pos ≈ SVector(30.5, 35.0)
    end

    @testset "No penetration" begin
        function worst_gap(D, rule; periodic, seed = 3, n = 20, steps = 200, R = 0.5)
            rng = Xoshiro(seed)
            L = 60.0
            space = ContinuousSpace(ntuple(_ -> L, D); periodic)
            spheres = Obstacles(space, [ObstacleSphere(rand(rng, SVector{D}) * L, 3 + 4rand(rng)) for _ in 1:15]; contact = rule)
            caps = Obstacles(space, [ObstacleCapsule(rand(rng, SVector{D}) * L, rand(rng, SVector{D}) * L, 0.3) for _ in 1:3]; contact = rule)
            boxes = Obstacles(space, [ObstacleBox(c, c .+ 4.0) for c in (rand(rng, SVector{D}) * (L - 4) for _ in 1:3)]; contact = rule)
            sets = periodic ? (spheres, caps, boxes) : (spheres, caps, boxes, Walls(space; contact = rule))
            model = StandardABM(Microbe{D,2}, space, 0.1; rng, properties = Dict(:obstacles => sets))
            motility = RunTumble(run_speed = [30.0], run_duration = 1.0, angle = Isotropic(D))
            for _ in 1:n
                add_agent!(random_free_position(model; radius = R), model; motility, radius = R, rotational_diffusivity = 0.1)
            end
            worst = Inf
            for _ in 1:steps
                run!(model, 1)
                for m in allagents(model)
                    worst = min(worst, MA._nearest(sets, m.pos, Inf)[1] - R)
                end
            end
            worst
        end
        for D in (2, 3), rule in (Halt(), Slide(), Slide(align = true), Reflect(), Randomize()), periodic in (false, true)
            @test worst_gap(D, rule; periodic) > -1e-9
        end

        # non-convex SDF (union of two disks): overlaps removed at the end of each step
        f(x) = min(norm(x - SVector(20.0, 30.0)) - 6, norm(x - SVector(28.0, 30.0)) - 6)
        space = ContinuousSpace((60.0, 60.0); periodic = false)
        sets = (Obstacles(space, [ObstacleSDF(f, (24.0, 30.0), 11.0)]; contact = Slide()), Walls(space))
        model = StandardABM(Microbe{2,2}, space, 0.1; rng = Xoshiro(1), properties = Dict(:obstacles => sets))
        motility = RunTumble(run_speed = [30.0], run_duration = 1.0, angle = Isotropic(2))
        for _ in 1:50
            add_agent!(random_free_position(model; radius = 0.5), model; motility, radius = 0.5, rotational_diffusivity = 0.1)
        end
        worst = Inf
        for _ in 1:300
            run!(model, 1)
            for m in allagents(model)
                worst = min(worst, MA._nearest(sets, m.pos, Inf)[1] - 0.5)
            end
        end
        @test worst > -1e-9
    end
end
