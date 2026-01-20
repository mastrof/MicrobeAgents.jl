using MicrobeAgents, Test
using LinearAlgebra: norm
using Random

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
        model = StandardABM(BrownBerg{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        add_agent!(model; motility, gain=600, receptor_binding_constant=100, memory=1)
        add_agent!(model; motility, gain=100, receptor_binding_constant=100, memory=1)
        add_agent!(model; motility, gain=600, receptor_binding_constant=10, memory=1)
        add_agent!(model; motility, gain=600, receptor_binding_constant=100, memory=2)
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1]) == bias(model[2]) == bias(model[3]) == bias(model[4]) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(BrownBerg{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, gain=600, receptor_binding_constant=100, memory=1)
        add_agent!(pos, model; vel, motility, gain=100, receptor_binding_constant=100, memory=1)
        add_agent!(pos, model; vel, motility, gain=600, receptor_binding_constant=5, memory=1)
        add_agent!(pos, model; vel, motility, gain=600, receptor_binding_constant=100, memory=2)
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1]) < bias(model[2])
        # maximum response occurs at K~C (5μM in this case)
        # so tumblebias smaller for the low K bacterium
        @test bias(model[1]) > bias(model[3])
        # longer memory => less affected by measurement => larger tumble bias
        @test bias(model[1]) < bias(model[4])

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = time_impulse_concentration,
            concentration_ramp = time_impulse_derivative,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(BrownBerg{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(model; motility, memory=1)
        add_agent!(model; motility, memory=2)
        adf, = run!(model, 5; adata=[bias])
        B = Analysis.adf_to_vectors(adf, :bias)
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
        @test bias(model[1]) ≈ 1
        @test bias(model[2]) ≈ 1
    end

    @testset "Brumley" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Brumley{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        add_agent!(model; motility, gain=5, memory=1, chemotactic_precision=0)
        add_agent!(model; motility, gain=1, memory=1, chemotactic_precision=0)
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1]) == bias(model[2]) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Brumley{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, gain=0.2, memory=1, chemotactic_precision=0)
        add_agent!(pos, model; vel, motility, gain=0.1, memory=1, chemotactic_precision=0)
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1]) < bias(model[2])
    end


    @testset "Celani" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Celani{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        add_agent!(model; motility, gain=5, memory=1)
        add_agent!(model; motility, gain=1, memory=1)
        add_agent!(model; motility, gain=5, memory=2)
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1]) == bias(model[2]) == bias(model[3]) == 1

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(Celani{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, gain=5, memory=1)
        add_agent!(pos, model; vel, motility, gain=1, memory=1)
        add_agent!(pos, model; vel, motility, gain=5, memory=2)
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1]) < bias(model[2])
        # longer memory => less affected by measurement => larger tumble bias
        @test bias(model[1]) < bias(model[3])
    end

    @testset "SonMenolascina" begin
        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = constant_background_concentration,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(SonMenolascina{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        add_agent!(model; motility, gain=600, memory=1)
        add_agent!(model; motility, gain=100, memory=1)
        add_agent!(model; motility, gain=600, memory=2)
        run!(model, 1)
        # no gradient => everyone has same bias and state independent of parameters
        @test bias(model[1]) == bias(model[2]) == bias(model[3]) == 1
        # since cT = 0.5 μM, speed should be 30% larger than specified
        @test speed(model[1]) == 20*1.3

        L = 100
        space = ContinuousSpace((L, L); periodic=false)
        dt = 0.1
        chemo = GenericChemoattractant{2}(;
            concentration_field = linear_x_concentration,
            concentration_gradient = linear_x_gradient,
        )
        properties = Dict(:chemoattractant => chemo)
        model = StandardABM(SonMenolascina{2,2}, space, dt; properties)
        motility = RunTumble(Inf, [20.0])
        pos = spacesize(model) ./ 2 # initialize at the center of domain
        vel = SVector(1.0, 0.0) # align on gradient direction
        add_agent!(pos, model; vel, motility, gain=600, memory=1)
        add_agent!(pos, model; vel, motility, gain=100, memory=1)
        add_agent!(pos, model; vel, motility, gain=600, memory=2)
        run!(model, 1)
        # larger gain => stronger response => longer runs => smaller tumble bias
        @test bias(model[1]) < bias(model[2])
        # longer memory => less affected by measurement => larger tumble bias
        @test bias(model[1]) < bias(model[3])
    end

end
