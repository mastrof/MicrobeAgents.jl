using MicrobeAgents, Test
using LinearAlgebra: norm
using Random

# adata columns are named after the accessor, so it must be a top-level function
const TUMBLE_MODEL = Ref{Any}()
tumblebias(m) = bias(m, TUMBLE_MODEL[])

const FLAT = (concentration_field = (m, model) -> 2.0,)
const XGRAD = (
    concentration_field = (m, model) -> position(m)[1] / 10,
    concentration_gradient = (m, model) -> SVector(1/10, 0.0),
)
chemo_model(field; dt = 0.1) = StandardABM(Microbe{2}, ContinuousSpace((100.0, 100.0); periodic = false), dt;
    properties = Dict(:chemoattractant => ChemicalField{2}(; field...)))
# microbe at the center of the domain, running up the x gradient
add_aligned!(model, behaviors; motility = RunTumble([20.0], Inf, Isotropic(2))) =
    add_agent!(SVector(50.0, 50.0), model; vel = SVector(1.0, 0.0), motility, behaviors)
# response of the chemotactic behavior alone (SonMenolascina also biases by speed)
chemobias(m, model) = bias(m.behaviors.chemotaxis, m, model)

@testset "Chemotaxis" verbose=true begin
    # (stronger, weaker) parameter pairs for each model
    gains = (
        BrownBerg = (BrownBerg(gain = 600, memory = 1), BrownBerg(gain = 100, memory = 1)),
        Brumley = (Brumley(gain = 0.2, memory = 1, chemotactic_precision = 0), Brumley(gain = 0.1, memory = 1, chemotactic_precision = 0)),
        Celani = (Celani(gain = 5, memory = 1), Celani(gain = 1, memory = 1)),
        SonMenolascina = (SonMenolascina(gain = 600, memory = 1), SonMenolascina(gain = 100, memory = 1)),
    )
    memories = (
        BrownBerg = (BrownBerg(gain = 600, memory = 1), BrownBerg(gain = 600, memory = 2)),
        Celani = (Celani(gain = 5, memory = 1), Celani(gain = 5, memory = 2)),
        SonMenolascina = (SonMenolascina(gain = 600, memory = 1), SonMenolascina(gain = 600, memory = 2)),
    )
    wrap(b) = b isa NamedTuple ? b : (chemotaxis = b,)

    @testset "$name: no response without gradients" for (name, pair) in pairs(gains)
        model = chemo_model(FLAT)
        foreach(b -> add_aligned!(model, wrap(b)), pair)
        run!(model, 1)
        @test chemobias(model[1], model) == chemobias(model[2], model) == 1
    end

    @testset "$name: larger gain, stronger response" for (name, (strong, weak)) in pairs(gains)
        model = chemo_model(XGRAD)
        add_aligned!(model, wrap(strong))
        add_aligned!(model, wrap(weak))
        run!(model, 1)
        # stronger response => longer runs up the gradient => smaller bias
        @test chemobias(model[1], model) < chemobias(model[2], model) < 1
    end

    @testset "$name: longer memory, weaker response" for (name, (short, long)) in pairs(memories)
        model = chemo_model(XGRAD)
        add_aligned!(model, wrap(short))
        add_aligned!(model, wrap(long))
        run!(model, 1)
        @test chemobias(model[1], model) < chemobias(model[2], model)
    end

    @testset "BrownBerg" begin
        # maximum response occurs at K ~ C (5μM in this case)
        model = chemo_model(XGRAD)
        add_aligned!(model, (chemotaxis = BrownBerg(gain = 600, receptor_binding_constant = 100),))
        add_aligned!(model, (chemotaxis = BrownBerg(gain = 600, receptor_binding_constant = 5),))
        run!(model, 1)
        @test chemobias(model[1], model) > chemobias(model[2], model)

        # response to a square pulse in time
        dt = 0.1
        model = chemo_model((
            concentration_field = (m, model) -> (1 <= abmtime(model) <= 3) ? 1.0 : 0.0,
            # positive dirac delta at t=1 and negative at t=3
            concentration_ramp = (m, model) -> abmtime(model) == 1 ? 1/dt : abmtime(model) == 3 ? -1/dt : 0.0,
        ); dt)
        add_aligned!(model, (BrownBerg(memory = 1),))
        add_aligned!(model, (BrownBerg(memory = 2),))
        TUMBLE_MODEL[] = model
        adf, = run!(model, 5; adata = [tumblebias])
        B = Analysis.adf_to_vectors(adf, :tumblebias)
        # positive response (bias < 1) when C goes up, negative (bias > 1) when C goes down;
        # longer memory gives a weaker response. Index 1 is timestep 0, hence the `1+`
        for b in B
            @test b[1+2] < 1 && b[1+3] < 1
            @test b[1+4] > 1 && b[1+5] > 1
        end
        @test abs.(B[1][1+2:end] .- 1) > abs.(B[2][1+2:end] .- 1)
        # after long time both return to null bias
        run!(model, 1000)
        @test bias(model[1], model) ≈ 1
        @test bias(model[2], model) ≈ 1
    end

    @testset "Noisy sensing requires a finite radius" begin
        model = chemo_model(FLAT)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = (Brumley(chemotactic_precision = 0),))
        run!(model, 5)
        @test bias(model[1], model) == 1
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = (Brumley(chemotactic_precision = 6),))
        m2 = add_agent!(model; motility, behaviors = (Celani(), Xie()))
        run!(model, 5)
        @test !isnan(bias(m2, model))
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = (Xie(chemotactic_precision = 1),))
        m3 = add_agent!(model; motility, radius = 0.5, behaviors = (Brumley(chemotactic_precision = 6),))
        run!(model, 5)
        @test isfinite(bias(m3, model))
    end

    @testset "Xie forward/backward gains" begin
        model = chemo_model(XGRAD)
        x = Xie(gain_forward = 2.7, gain_backward = 1.6)
        m = add_aligned!(model, (x,); motility = RunReverseFlick([20.0], Inf, [20.0], Inf))
        run!(model, 1)
        bf = bias(m, model) # forward run
        update_motilestate!(motilepattern(m), 3)
        bb = bias(m, model) # backward run
        @test bf != 1 && (bf - 1) / (bb - 1) ≈ 2.7 / 1.6
    end

    @testset "Chemokinesis" begin
        model = chemo_model(XGRAD) # C = x/10
        motility = RunTumble([20.0], Inf, Isotropic(2))
        low = add_agent!(SVector(5.0, 50.0), model; motility, vel = SVector(0.0, 1.0),
            behaviors = (Chemokinesis(threshold = 1.0, factor = 1.5),))
        high = add_agent!(SVector(50.0, 50.0), model; motility, vel = SVector(0.0, 1.0),
            behaviors = (Chemokinesis(threshold = 1.0, factor = 1.5),))
        run!(model, 1)
        @test speed(low) == 20.0
        @test speed(high) == 30.0
        @test low.speed == high.speed == 20.0 # base speed untouched
    end

    @testset "SpeedDependentTurnRate" begin
        model = chemo_model(FLAT)
        f(v) = 1 / (-0.55 / (1 + exp(-0.35 * (v - 18.88))) + 1.0)
        for v in (0.0, 10.0, 40.0)
            m = add_agent!(model; motility = RunTumble([v], 1.0, Isotropic(2)),
                behaviors = (SpeedDependentTurnRate(),))
            @test bias(m, model) ≈ f(v) / f(0)
        end
    end

    @testset "SonMenolascina" begin
        bs = SonMenolascina()
        @test keys(bs) == (:chemokinesis, :chemotaxis, :turnrate, :flick)
        @test length((; bs..., extra = (m, model) -> nothing)) == 5 # extendable
        # above the default threshold (0.05 μM) the speed is 30% larger
        model = chemo_model(FLAT)
        m = add_aligned!(model, SonMenolascina())
        run!(model, 1)
        @test speed(m) == 20 * 1.3

        # flick hook only acts when leaving the backward run of a 4-state motility
        model = chemo_model(FLAT)
        m = add_agent!(model; motility = RunReverseFlick([30.0], 1.0, [30.0], 1.0),
            behaviors = (flick = SpeedDependentFlick(),))
        w = [0.0, 0.0, 0.0, 1.0]
        update_motilestate!(motilepattern(m), 3)
        transition_weights!(w, m.behaviors.flick, m, model)
        p = 0.055 + 0.72 / (1 + exp(-0.25 * (30.0 - 36.0)))
        @test w ≈ [0.0, 1 - p, 0.0, p]
        w = [0.0, 0.0, 0.0, 1.0]
        update_motilestate!(motilepattern(m), 1)
        transition_weights!(w, m.behaviors.flick, m, model)
        @test w == [0.0, 0.0, 0.0, 1.0]
    end
end
