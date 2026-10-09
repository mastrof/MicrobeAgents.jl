using MicrobeAgents, Test, Random

@testset "Field cache" begin
    calls = Ref(0)
    cfield(m, model) = (calls[] += 1; position(m)[1])
    chemo = ChemicalField{2}(; concentration_field = cfield)
    space = ContinuousSpace((100.0, 100.0))
    model = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemicalfield => chemo))
    @test :field_cache in keys(abmproperties(model))
    reader(m, model) = MicrobeAgents.concentration(m, model)
    motility = RunTumble([1.0], Inf, Isotropic(2))
    add_agent!(SVector(10.0, 50.0), model; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader, reader, reader))
    calls[] = 0
    run!(model, 1)
    @test calls[] == 1 # three readers, one evaluation

    # outside a step: always fresh
    calls[] = 0
    @test MicrobeAgents.concentration(model[1], model) == 11.0
    @test MicrobeAgents.concentration(model[1], model) == 11.0
    @test calls[] == 2

    # reordered subroutines: affect before move must not leak stale values into reorient
    seen = Float64[]
    probe = Behavior(bias = (m, model) -> (push!(seen, MicrobeAgents.concentration(m, model)); 1.0))
    step_affect_first!(m, model) = (affect_step!(m, model); move_step!(m, model); reorient_step!(m, model))
    model2 = StandardABM(Microbe{2}, space, 1.0;
        properties = Dict(:chemicalfield => chemo), agent_step! = step_affect_first!)
    add_agent!(SVector(10.0, 50.0), model2; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader, probe))
    run!(model2, 1)
    @test seen == [11.0] # post-move position, not the cached pre-move 10.0

    # a behavior that moves the microbe invalidates the cache for later readers
    seen2 = Float64[]
    mover(m, model) = (m.pos = SVector(60.0, 50.0); MicrobeAgents.invalidate_field_cache!(model))
    reader2(m, model) = push!(seen2, MicrobeAgents.concentration(m, model))
    model3 = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemicalfield => chemo))
    add_agent!(SVector(10.0, 50.0), model3; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader2, mover, reader2))
    run!(model3, 1)
    @test seen2 == [11.0, 60.0]

    # type stability
    m = model[1]
    @test (@inferred MicrobeAgents.concentration(m, model)) isa Float64
    @test (@inferred gradient(m, model)) isa SVector{2,Float64}

    @testset "multiple fields" begin
        calls = Ref(0)
        fa(m, model) = (calls[] += 1; position(m)[1])
        fb(m, model) = (calls[] += 1; -position(m)[1])
        A = ChemicalField{2}(; concentration_field = fa)
        B = ChemicalField{2}(;
            concentration_field = fb,
            concentration_gradient = (m, model) -> SVector(-1.0, 0.0),
            diffusivity = (m, model) -> 100.0,
        )
        properties = Dict(:repellent => B, :chemicalfield => A)
        model = StandardABM(Microbe{2}, space, 1.0; properties)
        @test MicrobeAgents.check_field(model, :chemicalfield) === nothing
        @test MicrobeAgents.check_field(model, :repellent) === nothing
        err = try MicrobeAgents.check_field(model, :nope) catch e e end
        @test err isa ArgumentError
        @test occursin("repellent", err.msg) && occursin("chemicalfield", err.msg)

        reader_a(m, model) = MicrobeAgents.concentration(m, model, :chemicalfield)
        reader_b(m, model) = MicrobeAgents.concentration(m, model, :repellent)
        add_agent!(SVector(10.0, 50.0), model; motility, vel = SVector(1.0, 0.0),
            behaviors = (reader_a, reader_b, reader_a, reader_b))
        @test_throws ArgumentError MicrobeAgents.concentration(model[1], model, :nope)
        calls[] = 0
        run!(model, 1)
        @test calls[] == 2 # one evaluation per field per step

        m = model[1]
        # outside a step: fresh values of the requested field
        @test MicrobeAgents.concentration(m, model) == position(m)[1]
        @test MicrobeAgents.concentration(m, model, :repellent) == -position(m)[1]
        @test gradient(m, model, :repellent) == SVector(-1.0, 0.0)
        @test gradient(m, model) == SVector(0.0, 0.0)
        @test diffusivity(m, model, :repellent) == 100.0
        @test diffusivity(m, model) == 608.0

        # invalidation resets every field slot
        seen = Float64[]
        mover(m, model) = (m.pos = SVector(60.0, 50.0); MicrobeAgents.invalidate_field_cache!(model))
        rd(m, model) = push!(seen, MicrobeAgents.concentration(m, model, :repellent))
        model2 = StandardABM(Microbe{2}, space, 1.0; properties)
        add_agent!(SVector(10.0, 50.0), model2; motility, vel = SVector(1.0, 0.0),
            behaviors = (rd, mover, rd))
        run!(model2, 1)
        @test seen == [-11.0, -60.0]

        # cached path is allocation-free and type-stable for every field
        MicrobeAgents.reset_field_cache!(model, m)
        MicrobeAgents.concentration(m, model, :repellent)
        allocs(m, model, k) = @allocated MicrobeAgents.concentration(m, model, k)
        allocs(m, model, :repellent)
        @test allocs(m, model, :repellent) == 0
        @test (@inferred MicrobeAgents.concentration(m, model, :repellent)) isa Float64
        @test (@inferred gradient(m, model, :repellent)) isa SVector{2,Float64}
    end

    @testset "invalid field properties" begin
        space = ContinuousSpace((100.0, 100.0))
        e1 = try StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemicalfield => nothing)) catch e e end
        @test e1 isa ArgumentError && occursin("chemicalfield", e1.msg)
        e2 = try StandardABM(Microbe{2}, space, 1.0; properties = Dict(:x => ChemicalField{3}())) catch e e end
        @test e2 isa ArgumentError && occursin(":x", e2.msg)
    end
end
