using MicrobeAgents, Test, Random
include(joinpath(@__DIR__, "fixtures", "reference.jl"))

# same scenarios as fixtures/generate_reference.jl, with the behaviors API
const REF_L = 1000.0
const REF_DT = 0.1
const REF_NSTEPS = 300
const REF_NAGENTS = 3
ref_field(m, model) = 1.0 + 0.01 * position(m)[1] + 0.001 * abmtime(model) * REF_DT
ref_grad(m::AbstractMicrobe{D}, model) where {D} = SVector{D}(i == 1 ? 0.01 : 0.0 for i in 1:D)
ref_ramp(m, model) = 0.001
const REF_CHEMO = ChemicalField{2}(;
    concentration_field = ref_field,
    concentration_gradient = ref_grad,
    concentration_ramp = ref_ramp,
)
ref_rt() = RunTumble([30.0], 1.0, Isotropic(2); tumble_duration = 0.1)
ref_rrf() = RunReverseFlick([46.5], 0.45, [46.5], 0.45)

# adata columns are named after the accessor, so these must be top-level functions
const REF_MODEL = Ref{Any}()
xpos(m) = position(m)[1]
ypos(m) = position(m)[2]
chemotaxis_bias(m) = bias(m.behaviors.chemotaxis, m, REF_MODEL[])
chemotaxis_state(m) = m.behaviors.chemotaxis.state
state_m(m) = m.behaviors.chemotaxis.state_m
state_z(m) = m.behaviors.chemotaxis.state_z

function run_scenario(motility, behaviors, kw, extras)
    space = ContinuousSpace((REF_L, REF_L); periodic = true)
    model = StandardABM(Microbe{2}, space, REF_DT;
        properties = Dict(:chemoattractant => REF_CHEMO),
        rng = Xoshiro(1234), container = Vector,
    )
    for _ in 1:REF_NAGENTS
        add_agent!(model; motility = motility(), behaviors, kw...)
    end
    REF_MODEL[] = model
    adata = [xpos, ypos, speed, chemotaxis_bias, chemotaxis_state, extras...]
    adf, = run!(model, REF_NSTEPS; adata, when = 0:10:REF_NSTEPS)
    return adf
end

function test_against_reference(name, adf)
    ref = REFERENCE[name]
    @testset "$name reproduces reference" begin
        for k in keys(ref)
            @test Vector{Float64}(adf[!, k]) == ref[k]
        end
    end
end

@testset "Regression vs pre-behaviors API" begin
    test_against_reference("BrownBerg", run_scenario(ref_rt,
        (chemotaxis = BrownBerg(),),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
    test_against_reference("Brumley", run_scenario(ref_rrf,
        (chemotaxis = Brumley(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
    test_against_reference("Celani", run_scenario(ref_rt,
        (chemotaxis = Celani(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.26, radius = 0.5), []))
    test_against_reference("Xie", run_scenario(ref_rrf,
        (chemotaxis = Xie(chemotactic_precision = 6.0),),
        (rotational_diffusivity = 0.26, radius = 0.5), [state_m, state_z]))
    test_against_reference("SonMenolascina", run_scenario(
        () -> RunReverseFlick([30.0], 0.5, [30.0], 0.5),
        SonMenolascina(),
        (rotational_diffusivity = 0.035, radius = 0.5), []))
end
