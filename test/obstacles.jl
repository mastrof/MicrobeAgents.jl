using MicrobeAgents, Test
using LinearAlgebra: norm, dot
using Random
using MicrobeAgents.Agents: allagents

# one microbe moving straight (no rotational diffusion, never tumbles)
function straight_model(space, obstacles, pos, vel; speed = 10.0, radius = 0.5, dt = 1.0)
    D = length(pos)
    model = StandardABM(Microbe{D}, space, dt; rng = Xoshiro(1),
        properties = Dict(:obstacles => obstacles))
    motility = RunTumble(run_speed = [speed], run_duration = Inf, angle = Isotropic(D))
    m = add_agent!(SVector{D}(pos), model; motility, vel = SVector{D}(vel), radius)
    return model, m
end

# closest surface from a point (microbe of zero radius)
function query(space, obstacles, x; cutoff = Inf)
    model, m = straight_model(space, obstacles, x, ntuple(i -> i == 1 ? 1.0 : 0.0, length(x)); radius = 0.0)
    return nearest_surface(m, model; cutoff)
end

@testset "Obstacles" begin
    @testset "Primitives" begin
        @test sdf(ObstacleSphere((0.0, 0.0), 2.0), SVector(3.0, 4.0)) ≈ 3.0
        @test sdf(ObstacleHalfSpace((0.0, 0.0), (0.0, 2.0)), SVector(5.0, 3.0)) ≈ 3.0
        b = ObstacleBox((0.0, 0.0), (2.0, 4.0))
        @test sdf(b, SVector(3.0, 2.0)) ≈ 1.0
        @test sdf(b, SVector(1.0, 2.0)) ≈ -1.0
        @test sdf(b, SVector(3.0, 6.0)) ≈ sqrt(5)
        @test_throws ArgumentError ObstacleBox((2.0, 0.0), (0.0, 1.0))
        @test sdf(ObstacleBox(SVector(1.0, 1.0), SVector(3.0, 5.0)), SVector(4.0, 3.0)) ≈ 1.0 # SVector args
        c = ObstacleCapsule((0.0, 0.0), (4.0, 0.0), 1.0)
        @test sdf(c, SVector(2.0, 3.0)) ≈ 2.0
        @test sdf(c, SVector(6.0, 0.0)) ≈ 1.0
        @test sdf(ObstacleSDF(x -> norm(x) - 1, (0.0, 0.0), 1.0), SVector(2.0, 0.0)) ≈ 1.0
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
        @test contact!(Slide(), m, c, model) ≈ SVector(0.0, 4.0)
        @test direction(m) == SVector(-0.6, 0.8)
        away = SurfaceContact{2}(1, 1, SVector(24.5, 25.0), n, SVector(3.0, 4.0))
        @test contact!(Slide(), m, away, model) ≈ SVector(3.0, 4.0) # not into the surface
        @test contact!(Slide(align = true), m, c, model) ≈ SVector(0.0, 4.0)
        @test direction(m) ≈ SVector(0.0, 1.0)

        m.vel = SVector(-0.6, 0.8)
        @test contact!(Reflect(), m, c, model) ≈ SVector(3.0, 4.0)
        @test direction(m) ≈ SVector(0.6, 0.8)

        m.vel = SVector(-0.6, 0.8)
        u = contact!(Randomize(), m, c, model)
        @test dot(direction(m), n) ≥ 0
        @test u ≈ norm(s) * direction(m)

        @test contact!((mic, con, mod) -> 2 .* con.displacement, m, c, model) == 2 .* s
    end

    @testset "Obstacle sets and nearest_surface" begin
        # compare with a brute-force signed distance over all primitives and periodic images
        rng = Xoshiro(1)
        for D in 2:3, periodic in (false, true), (cellsize, margin) in ((nothing, nothing), (3.0, 0.5))
            L = 50.0
            space = ContinuousSpace(ntuple(_ -> L, D); periodic)
            prims = [ObstacleSphere(rand(rng, SVector{D}) * L, 1 + 2rand(rng)) for _ in 1:30]
            set = Obstacles(space, prims; cellsize, margin)
            model, m = straight_model(space, (set,), fill(1.0, D), ntuple(i -> i == 1 ? 1.0 : 0.0, D); radius = 0.0)
            images = Iterators.product(ntuple(d -> periodic ? (-L, 0.0, L) : (0.0,), D)...)
            @test all(1:25) do _
                x = rand(rng, SVector{D}) * L
                m.pos = x
                cutoff = rand(rng, (0.3, 3.0, 12.0, Inf))
                b = minimum(sdf(p, x - SVector{D}(s)) for p in prims for s in images)
                q = nearest_surface(m, model; cutoff)
                # includes points inside solids
                b ≤ cutoff ? q.distance ≈ b : (q.distance == Inf && q.index == 0)
            end
        end

        # surface normals of the primitives
        space = ContinuousSpace((50.0, 50.0); periodic = false)
        for (p, x, n) in (
            (ObstacleSphere((20.0, 20.0), 3.0), (26.0, 28.0), (0.6, 0.8)),
            (ObstacleBox((20.0, 20.0), (22.0, 24.0)), (23.0, 26.0), (1.0, 2.0) ./ sqrt(5)),
            (ObstacleCapsule((20.0, 20.0), (24.0, 20.0), 1.0), (22.0, 23.0), (0.0, 1.0)),
            (ObstacleSDF(x -> norm(x - SVector(20.0, 20.0)) - 3, (20.0, 20.0), 3.0), (26.0, 28.0), (0.6, 0.8)),
        )
            @test query(space, (Obstacles(space, [p]),), x).normal ≈ SVector(n) atol = 1e-6
        end

        # obstacle straddling a periodic boundary is seen from the other side
        pspace = ContinuousSpace((50.0, 50.0); periodic = true)
        q = query(pspace, (Obstacles(pspace, [ObstacleSphere((49.0, 25.0), 3.0)]),), (5.0, 25.0))
        @test q.distance ≈ 3.0 && q.normal ≈ SVector(1.0, 0.0)

        # walls and tuples of sets
        space = ContinuousSpace((50.0, 40.0); periodic = (true, false))
        walls = Walls(space)
        sets = (Obstacles(space, [ObstacleSphere((25.0, 20.0), 5.0)]), walls)
        q = query(space, sets, (10.0, 3.0))
        @test q.distance ≈ 3.0 && q.normal ≈ SVector(0.0, 1.0) && q.set == 2 && q.index == 1
        q = query(space, sets, (10.0, 39.0))
        @test q.distance ≈ 1.0 && q.normal ≈ SVector(0.0, -1.0) && q.set == 2 && q.index == 2
        q = query(space, sets, (25.0, 12.0))
        @test q.distance ≈ 3.0 && q.set == 1
        @test query(space, sets, (25.0, 12.0); cutoff = 1.0).distance == Inf
        @test_throws ArgumentError Walls(ContinuousSpace((50.0, 50.0); periodic = true))

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
        @test abmproperties(StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => sets))).obstacles === sets
        other = ContinuousSpace((60.0, 50.0); periodic = false)
        @test_throws ArgumentError StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => (Walls(other),)))
        @test_throws ArgumentError StandardABM(Microbe{2}, space, 0.1;
            properties = Dict(:obstacles => [walls]))

        # nearest_surface accounts for the microbe radius; no obstacles -> Inf
        motility = RunTumble(run_speed = [10.0], run_duration = 1.0, angle = Isotropic(2))
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:obstacles => sets))
        m = add_agent!(SVector(8.0, 25.0), model; motility, radius = 0.5)
        @test nearest_surface(m, model).distance ≈ 7.5
        empty_model = StandardABM(Microbe{2}, space, 0.1)
        m = add_agent!(SVector(10.0, 25.0), empty_model; motility)
        @test nearest_surface(m, empty_model).distance == Inf

        # random_free_position
        model = StandardABM(Microbe{2}, space, 0.1; rng = Xoshiro(2), properties = Dict(:obstacles => sets))
        @test all(1:5) do _
            x = random_free_position(model; radius = 1.0)
            norm(x - SVector(25.0, 25.0)) ≥ 6.0 && all(1.0 .≤ x .≤ 49.0)
        end
        full = (Obstacles(space, [ObstacleBox((-1.0, -1.0), (51.0, 51.0))]),)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:obstacles => full))
        @test_throws ErrorException random_free_position(model; maxtries = 10)
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

        # obstacle across the periodic boundary, approached from the other side
        pspace = ContinuousSpace((50.0, 50.0); periodic = true)
        model, m = straight_model(pspace, (Obstacles(pspace, [ObstacleSphere((49.0, 25.0), 3.0)]),),
            (5.0, 25.0), (-1.0, 0.0))
        run!(model, 1)
        @test m.pos ≈ SVector(2.5, 25.0)

        # a microbe placed overlapping an obstacle is pushed out on its first step
        disk = (Obstacles(space, [ObstacleSphere((25.0, 25.0), 5.0)]),)
        model, m = straight_model(space, disk, (27.0, 25.0), (0.0, 1.0))
        run!(model, 1)
        @test m.pos ≈ SVector(30.5, 25.0)
        run!(model, 1) # then leaves along the tangent
        @test m.pos ≈ SVector(30.5, 35.0)

        # Slide makes progress along curved and oblique surfaces (round-off in the
        # normal component of the projected displacement must not count as contact)
        c0 = SVector(25.0, 25.0)
        for θ in (0.3, 4.1)
            n = SVector(cos(θ), sin(θ)); t = SVector(-sin(θ), cos(θ))
            d = (t - 0.3n) / norm(t - 0.3n)
            sets = (Obstacles(space, [ObstacleSphere(c0, 5.0)]; contact = Slide()),)
            model, m = straight_model(space, sets, c0 + 5.5n, d; dt = 0.1)
            x0 = m.pos
            run!(model, 1)
            @test norm(m.pos - x0) > 0.5
            sets = (Obstacles(space, [ObstacleHalfSpace(c0, -n)]; contact = Slide()),)
            model, m = straight_model(space, sets, c0 - 0.5n, -d; dt = 0.1)
            x0 = m.pos
            run!(model, 1)
            @test norm(m.pos - x0) > 0.5
        end

        # Slide(align = true) keeps following a curved surface until reoriented
        for (sp, c0, t0) in ((space, SVector(25.0, 25.0), SVector(0.0, 1.0)),
                (ContinuousSpace((50.0, 50.0); periodic = true), SVector(1.0, 25.0), SVector(0.0, -1.0)))
            n0 = SVector(1.0, 0.0)
            sets = (Obstacles(sp, [ObstacleSphere(c0, 5.0)]; contact = Slide(align = true)),)
            model, m = straight_model(sp, sets, c0 + 5.5n0, -0.3n0 + t0; dt = 0.1)
            @test all(1:10) do _
                run!(model, 1)
                abs(nearest_surface(m, model).distance) < 1e-6
            end
            m.vel = nearest_surface(m, model).normal # reoriented away from the surface
            run!(model, 5)
            @test nearest_surface(m, model).distance > 1.0
        end

        # unindexed (unbounded) obstacles see periodic images
        pspace100 = ContinuousSpace((100.0, 100.0); periodic = true)
        unb = (Obstacles(pspace100, [ObstacleSDF(x -> norm(x - SVector(1.0, 50.0)) - 5, (1.0, 50.0))]),)
        model, m = straight_model(pspace100, unb, (93.0, 50.0), (1.0, 0.0); dt = 0.1)
        @test nearest_surface(m, model).distance ≈ 2.5 atol = 1e-6
        run!(model, 6)
        @test m.pos ≈ SVector(95.5, 50.0) atol = 1e-6
    end

    @testset "No penetration" begin
        # many microbes bouncing among random spheres, capsules and boxes
        function worst_gap(D, rule, periodic; R = 0.5)
            rng = Xoshiro(3)
            L = 60.0
            space = ContinuousSpace(ntuple(_ -> L, D); periodic)
            spheres = Obstacles(space, [ObstacleSphere(rand(rng, SVector{D}) * L, 3 + 4rand(rng)) for _ in 1:15]; contact = rule)
            caps = Obstacles(space, [ObstacleCapsule(rand(rng, SVector{D}) * L, rand(rng, SVector{D}) * L, 0.3) for _ in 1:3]; contact = rule)
            boxes = Obstacles(space, [ObstacleBox(c, c .+ 4.0) for c in (rand(rng, SVector{D}) * (L - 4) for _ in 1:3)]; contact = rule)
            sets = periodic ? (spheres, caps, boxes) : (spheres, caps, boxes, Walls(space; contact = rule))
            model = StandardABM(Microbe{D,2}, space, 0.1; rng, properties = Dict(:obstacles => sets))
            motility = RunTumble(run_speed = [30.0], run_duration = 1.0, angle = Isotropic(D))
            for _ in 1:10
                add_agent!(random_free_position(model; radius = R), model; motility, radius = R, rotational_diffusivity = 0.1)
            end
            worst = Inf
            for _ in 1:60
                run!(model, 1)
                for m in allagents(model)
                    worst = min(worst, nearest_surface(m, model).distance)
                end
            end
            worst
        end
        for (D, rule, periodic) in ((2, Halt(), false), (2, Slide(), true), (2, Slide(align = true), false),
                (3, Reflect(), false), (3, Randomize(), true))
            @test worst_gap(D, rule, periodic) > -1e-9
        end

        # non-convex SDF (union of two disks): overlaps removed at the end of each step
        f(x) = min(norm(x - SVector(20.0, 30.0)) - 6, norm(x - SVector(28.0, 30.0)) - 6)
        space = ContinuousSpace((60.0, 60.0); periodic = false)
        sets = (Obstacles(space, [ObstacleSDF(f, (24.0, 30.0), 11.0)]; contact = Slide()), Walls(space))
        model = StandardABM(Microbe{2,2}, space, 0.1; rng = Xoshiro(1), properties = Dict(:obstacles => sets))
        motility = RunTumble(run_speed = [30.0], run_duration = 1.0, angle = Isotropic(2))
        for _ in 1:15
            add_agent!(random_free_position(model; radius = 0.5), model; motility, radius = 0.5, rotational_diffusivity = 0.1)
        end
        worst = Inf
        for _ in 1:80
            run!(model, 1)
            for m in allagents(model)
                worst = min(worst, nearest_surface(m, model).distance)
            end
        end
        @test worst > -1e-9
    end

    @testset "Non-convex tracing" begin
        # grazing incidence on a flat wall: a non-convex SDF must match the half-space
        space = ContinuousSpace((200.0, 100.0); periodic = false)
        plane = ObstacleHalfSpace((0.0, 50.0), (0.0, 1.0))
        wall = ObstacleSDF(x -> x[2] - 50, (100.0, 25.0))
        for θ in (2, 10)
            v = (cosd(θ), -sind(θ))
            model1, m1 = straight_model(space, (Obstacles(space, [plane]; contact = Reflect()),), (10.0, 51.5), v; speed = 30.0, dt = 0.1)
            model2, m2 = straight_model(space, (Obstacles(space, [wall]; contact = Reflect()),), (10.0, 51.5), v; speed = 30.0, dt = 0.1)
            @test all(1:40) do _
                run!(model1, 1)
                run!(model2, 1)
                norm(m2.pos - m1.pos) < 1e-4
            end
        end

        # sliding inside a thin ring must not tunnel through it
        space = ContinuousSpace((100.0, 100.0); periodic = false)
        c = SVector(50.0, 50.0)
        ring = ObstacleSDF(x -> abs(norm(x - c) - 10) - 0.1, c, 10.1)
        for rule in (Slide(), Slide(align = true))
            model, m = straight_model(space, (Obstacles(space, [ring]; contact = rule),),
                c + SVector(0.0, 9.4), (1.0, 0.0); speed = 60.0, dt = 0.1)
            @test all(1:50) do _
                run!(model, 1)
                norm(m.pos - c) < 9.9 - 0.5 + 1e-6
            end
        end
    end
end
