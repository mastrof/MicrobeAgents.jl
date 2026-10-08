using MicrobeAgents, Test
using LinearAlgebra: norm
using Random

# adata columns are named after the accessor, so it must be a top-level function
const TUMBLE_MODEL = Ref{Any}()
tumblebias(m) = bias(m, TUMBLE_MODEL[])

@testset "Chemotaxis" verbose=true begin
    function constant_background_concentration(microbe, model)
        2.0
    end
    function linear_x_concentration(microbe, model)
        position(microbe)[1] / 10
    end
    function linear_x_gradient(microbe::AbstractMicrobe{D,N}, model) where {D,N}
        SVector{D}(i == 1 ? 1/10 : 0.0 for i in 1:D)
    end
    function time_impulse_concentration(microbe, model)
        t = abmtime(model)
        # square pulse for 3 timesteps
        (1 <= t <= 3) ? 1.0 : 0.0
    end
    function time_impulse_derivative(microbe, model)
        t = abmtime(model)
        dt = model.timestep
        # positive dirac delta at t=1 and negative at t=3
        if t == 1
            1.0 / dt
        elseif t == 3
            -1.0 /dt
        else
            0.0
        end

    end

    @testset "BrownBerg" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=100, memory=1),))
        add_agent!(model; motility, behaviors = (BrownBerg(gain=100, receptor_binding_constant=100, memory=1),))
        add_agent!(model; motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=10, memory=1),))
        add_agent!(model; motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=100, memory=2),))
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1], model) == bias(model[2], model) == bias(model[3], model) == bias(model[4], model) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=100, memory=1),))
        add_agent!(pos, model; vel, motility, behaviors = (BrownBerg(gain=100, receptor_binding_constant=100, memory=1),))
        add_agent!(pos, model; vel, motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=5, memory=1),))
        add_agent!(pos, model; vel, motility, behaviors = (BrownBerg(gain=600, receptor_binding_constant=100, memory=2),))
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1], model) < bias(model[2], model)
        # maximum response occurs at K~C (5μM in this case)
        # so tumblebias smaller for the low K bacterium
        @test bias(model[1], model) > bias(model[3], model)
        # longer memory => less affected by measurement => larger tumble bias
        @test bias(model[1], model) < bias(model[4], model)

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = time_impulse_concentration,
            concentration_ramp = time_impulse_derivative,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(model; motility, behaviors = (BrownBerg(memory=1),))
        add_agent!(model; motility, behaviors = (BrownBerg(memory=2),))
        TUMBLE_MODEL[] = model
        adf, = run!(model, 5; adata=[tumblebias])
        B = Analysis.adf_to_vectors(adf, :tumblebias)
        # positive response (bias < 1) when C goes up
        # negative response (bias > 1) when C goes down
        # bacterium with longer memory has weaker response (closer to 1)
        # time index 1 corresponds to timestep 0 hence the `1+`
        @test (B[1][1+2] < 1) && (B[1][1+3] < 1)
        @test B[2][1+2] < 1 && B[2][1+3] < 1
        @test B[1][1+4] > 1 && B[1][1+5] > 1
        @test B[2][1+4] > 1 && B[2][1+5] > 1
        @test abs.(B[1][1+2:end] .- 1) > abs.(B[2][1+2:end] .- 1)
        # after long time both of them return to null bias
        run!(model, 1000)
        @test bias(model[1], model) ≈ 1
        @test bias(model[2], model) ≈ 1
    end

    @testset "Brumley" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = (Brumley(gain=5, memory=1, chemotactic_precision=0),))
        add_agent!(model; motility, behaviors = (Brumley(gain=1, memory=1, chemotactic_precision=0),))
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1], model) == bias(model[2], model) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, behaviors = (Brumley(gain=0.2, memory=1, chemotactic_precision=0),))
        add_agent!(pos, model; vel, motility, behaviors = (Brumley(gain=0.1, memory=1, chemotactic_precision=0),))
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1], model) < bias(model[2], model)
    end


    @testset "Noise-free sensing with zero radius" begin
        space = ContinuousSpace((100.0, 100.0))
        chemo = GenericChemoattractant{2}(; concentration_field = constant_background_concentration)
        model = StandardABM(Microbe{2}, space, 0.1; properties = Dict(:chemoattractant => chemo))
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
    end

    @testset "Xie forward/backward gains" begin
        space = ContinuousSpace((100.0, 100.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        motility = RunReverseFlick([0.0], Inf, [0.0], Inf)
        add_agent!(model; motility, behaviors = (Xie(),))
        m = model[1]
        x = m.behaviors[1]
        x.state = 0.5
        @test bias(m, model) == 1 + x.gain_forward * 0.5
        m.motility.current_state = 3
        @test bias(m, model) == 1 + x.gain_backward * 0.5
    end

    @testset "Celani" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = (Celani(gain=5, memory=1),))
        add_agent!(model; motility, behaviors = (Celani(gain=1, memory=1),))
        add_agent!(model; motility, behaviors = (Celani(gain=5, memory=2),))
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1], model) == bias(model[2], model) == bias(model[3], model) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, behaviors = (Celani(gain=5, memory=1),))
        add_agent!(pos, model; vel, motility, behaviors = (Celani(gain=1, memory=1),))
        add_agent!(pos, model; vel, motility, behaviors = (Celani(gain=5, memory=2),))
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1], model) < bias(model[2], model)
        # longer memory => less affected by measurement => larger tumble bias
        @test bias(model[1], model) < bias(model[3], model)
    end

    @testset "SonMenolascina" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        motility = RunTumble([20.0], Inf, Isotropic(2))
        add_agent!(model; motility, behaviors = SonMenolascina(gain=600, memory=1))
        add_agent!(model; motility, behaviors = SonMenolascina(gain=100, memory=1))
        add_agent!(model; motility, behaviors = SonMenolascina(gain=600, memory=2))
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1].behaviors.chemotaxis, model[1], model) ==
              bias(model[2].behaviors.chemotaxis, model[2], model) ==
              bias(model[3].behaviors.chemotaxis, model[3], model) == 1
        # since threshold = 0.05 μM, speed should be 30% larger than specified
        @test speed(model[1]) == 20*1.3

        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Microbe{2}, space, dt; properties)
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, behaviors = SonMenolascina(gain=600, memory=1))
        add_agent!(pos, model; vel, motility, behaviors = SonMenolascina(gain=100, memory=1))
        add_agent!(pos, model; vel, motility, behaviors = SonMenolascina(gain=600, memory=2))
        run!(model, 1)
        b(i) = bias(model[i].behaviors.chemotaxis, model[i], model)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test b(1) < b(2)
        # longer memory => less affected by measurement => larger tumble bias
        @test b(1) < b(3)
    end

    @testset "Composable SonMenolascina pieces" begin
        bs = SonMenolascina()
        @test keys(bs) == (:chemokinesis, :chemotaxis, :turnrate, :flick)
        extended = (; bs..., extra = (m, model) -> nothing)
        @test length(extended) == 5
        # flick hook only acts when leaving the backward run of a 4-state motility
        space = ContinuousSpace((100.0, 100.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        add_agent!(model; motility = RunReverseFlick([30.0], 1.0, [30.0], 1.0),
            behaviors = (flick = SpeedDependentFlick(),))
        m = model[1]
        w = [0.0, 0.0, 0.0, 1.0]
        m.motility.current_state = 3
        transition_weights!(w, m.behaviors.flick, m, model)
        p = 0.055 + 0.72 / (1 + exp(-0.25 * (30.0 - 36.0)))
        @test w ≈ [0.0, 1 - p, 0.0, p]
        w = [0.0, 0.0, 0.0, 1.0]
        m.motility.current_state = 1
        transition_weights!(w, m.behaviors.flick, m, model)
        @test w == [0.0, 0.0, 0.0, 1.0]
    end
end
