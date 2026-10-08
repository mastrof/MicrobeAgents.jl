using MicrobeAgents, Test, Random

@testset "Field cache" begin
    calls = Ref(0)
    cfield(m, model) = (calls[] += 1; position(m)[1])
    chemo = GenericChemoattractant{2}(; concentration_field = cfield)
    space = ContinuousSpace((100.0, 100.0))
    model = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemoattractant => chemo))
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
        properties = Dict(:chemoattractant => chemo), agent_step! = step_affect_first!)
    add_agent!(SVector(10.0, 50.0), model2; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader, probe))
    run!(model2, 1)
    @test seen == [11.0] # post-move position, not the cached pre-move 10.0

    # a behavior that moves the microbe invalidates the cache for later readers
    seen2 = Float64[]
    mover(m, model) = (m.pos = SVector(60.0, 50.0); MicrobeAgents.invalidate_field_cache!(model))
    reader2(m, model) = push!(seen2, MicrobeAgents.concentration(m, model))
    model3 = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemoattractant => chemo))
    add_agent!(SVector(10.0, 50.0), model3; motility, vel = SVector(1.0, 0.0),
        behaviors = (reader2, mover, reader2))
    run!(model3, 1)
    @test seen2 == [11.0, 60.0]

    # type stability
    m = model[1]
    @test (@inferred MicrobeAgents.concentration(m, model)) isa Float64
    @test (@inferred gradient(m, model)) isa SVector{2,Float64}
end
