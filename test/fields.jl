using MicrobeAgents, Test, Random

# custom chemical field: type and methods must be defined at top level
struct XRamp <: AbstractChemicalField{2}
    slope::Float64
end
MicrobeAgents.concentration(c::XRamp) = (m, model) -> c.slope * position(m)[1]
MicrobeAgents.gradient(c::XRamp) = (m, model) -> SVector(c.slope, 0.0)
MicrobeAgents.time_derivative(c::XRamp) = (m, model) -> 0.0
MicrobeAgents.chemoattractant_diffusivity(c::XRamp) = (m, model) -> 500.0

@testset "Chemical fields" begin
    space = ContinuousSpace((100.0, 100.0))
    motility = RunTumble([1.0], Inf, Isotropic(2))

    @testset "Defaults" begin
        model = StandardABM(Microbe{2}, space, 1.0)
        m = add_agent!(SVector(10.0, 50.0), model; motility)
        @test chemoattractant(model) isa ChemicalField{2}
        @test concentration(m, model) === 0.0
        @test gradient(m, model) === SVector(0.0, 0.0)
        @test time_derivative(m, model) === 0.0
        @test chemoattractant_diffusivity(m, model) === 608.0
    end

    @testset "ChemicalField accessors" begin
        chemo = ChemicalField{2}(;
            concentration_field = (m, model) -> position(m)[1],
            concentration_gradient = (m, model) -> SVector(1.0, 0.0),
            concentration_ramp = (m, model) -> 0.5,
            diffusivity = (m, model) -> 100.0,
        )
        model = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemoattractant => chemo))
        m = add_agent!(SVector(10.0, 50.0), model; motility)
        @test chemoattractant(model) === chemo
        @test concentration(m, model) == 10.0
        @test concentration(model)(m, model) == 10.0
        @test gradient(m, model) == SVector(1.0, 0.0)
        @test time_derivative(m, model) == 0.5
        @test chemoattractant_diffusivity(m, model) == 100.0
        @test (@inferred concentration(m, model)) isa Float64
        @test (@inferred gradient(m, model)) isa SVector{2,Float64}
        # field functions must return Float64 / SVector{D,Float64}
        bad = ChemicalField{2}(; concentration_field = (m, model) -> 1)
        model = StandardABM(Microbe{2}, space, 1.0; properties = Dict(:chemoattractant => bad))
        m = add_agent!(model; motility)
        @test_throws TypeError concentration(m, model)
    end

    @testset "Custom AbstractChemicalField" begin
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:chemoattractant => XRamp(0.2)))
        vel = SVector(1.0, 0.0)
        up = add_agent!(SVector(50.0, 50.0), model; motility, vel, behaviors = (BrownBerg(),))
        down = add_agent!(SVector(50.0, 50.0), model; motility, vel = -vel, behaviors = (BrownBerg(),))
        @test concentration(up, model) == 10.0
        @test gradient(up, model) == SVector(0.2, 0.0)
        @test time_derivative(up, model) == 0.0
        @test chemoattractant_diffusivity(up, model) == 500.0
        run!(model, 1)
        # the behavior reads the custom field: up-gradient runs are extended
        @test bias(up, model) < 1 < bias(down, model)
    end

    @testset "Values within a step are fresh" begin
        chemo = ChemicalField{2}(; concentration_field = (m, model) -> position(m)[1])
        properties = Dict(:chemoattractant => chemo)
        vel = SVector(1.0, 0.0)
        # reordered subroutines: affect before move must not leak stale values into reorient
        seen = Float64[]
        probe = Behavior(bias = (m, model) -> (push!(seen, concentration(m, model)); 1.0))
        reader(m, model) = concentration(m, model)
        step_affect_first!(m, model) = (affect_step!(m, model); move_step!(m, model); reorient_step!(m, model))
        model = StandardABM(Microbe{2}, space, 1.0; properties, agent_step! = step_affect_first!)
        add_agent!(SVector(10.0, 50.0), model; motility, vel, behaviors = (reader, probe))
        run!(model, 1)
        @test seen == [11.0] # post-move position, not the pre-move 10.0

        # a behavior that moves the microbe invalidates the cache for later readers
        seen = Float64[]
        mover(m, model) = (m.pos = SVector(60.0, 50.0); MicrobeAgents.invalidate_field_cache!(model))
        reader2(m, model) = push!(seen, concentration(m, model))
        model = StandardABM(Microbe{2}, space, 1.0; properties)
        add_agent!(SVector(10.0, 50.0), model; motility, vel, behaviors = (reader2, mover, reader2))
        run!(model, 1)
        @test seen == [11.0, 60.0]
    end
end
