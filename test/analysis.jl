using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

@testset "Post-processing and analysis routines" begin
    dt = 1
    L = 100
    for D in 1:3
        extent = fill(float(L), SVector{D})
        space = ContinuousSpace(extent)
        pos = zero(SVector{D})
        U = 1
        #motility = RunTumble(speed=[U])
        turn_rate = 0
        motility = RunTumble(1/turn_rate, [U], 0.0)
        vel = fill(1/√D, SVector{D})
        model = StandardABM(Microbe{D}, space, dt)
        add_agent!(pos, model; vel, motility)
        nsteps = 10
        adata = [position]
        adf, mdf = run!(model, nsteps; adata)
        traj = Analysis.adf_to_matrix(adf, :position)
        @test traj isa AbstractMatrix{<:SVector}
        @test size(traj) == (nsteps+1, 1)
        real_traj = [pos .+ vel.*(U*n*dt) for n in 0:nsteps, _ in 1:1]
        Δ = norm.([traj[i] .- real_traj[i] for i in eachindex(traj)])
        for i in eachindex(Δ)
            @test Δ[i] ≈ 0 atol=1e-12
        end
    end

    #=
    @testset "Turn detection and run statistics" begin
        dt = 0.1
        L = 100
        for D in 1:3
            extent = fill(float(L), SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt)
            add_agent!(model; turn_rate=Inf) # turns every step
            nsteps = 10
            adf, _ = run!(model, nsteps; adata=[velocity])
            Analysis.detect_turns!(adf) # by default the new key is :has_turned
            println(names(adf))
            @test hasproperty(adf, :has_turned)
            @test Analysis.run_durations(adf) == [repeat([1],nsteps-1)]

            D == 1 && continue # only test rotational diffusion for D>1
            model = StandardABM(Microbe{D}, space, dt)
            # never turns but makes small deviations due to brownian noise
            Drot = 0.1
            add_agent!(model; turn_rate=0.0, rotational_diffusivity=Drot)
            nsteps = 10
            adf, _ = run!(model, nsteps; adata=[:vel])
            # detects turns at each step if threshold is not set
            Analysis.detect_turns!(adf)
            @test adf.has_turned == ones(Bool, length(adf.has_turned))
            # set threshold at 4σ to ignore rotational diffusion
            threshold_angle = 4 * sqrt(2*Drot*dt)
            Analysis.detect_turns!(adf; threshold_angle)
            @test adf.has_turned == zeros(Bool, length(adf.has_turned))
            @test Analysis.run_durations(adf) == [Float64[]]
        end
    end
    =#

    @testset "Autocorrelation function" begin
        for D in 1:3
            dt = 0.1
            nsteps = 100
            adata = [direction]
            L = 100
            extent = fill(float(L), SVector{D})
            space = ContinuousSpace(extent)
            rng = Xoshiro(35)
            model_periodic = StandardABM(Microbe{D}, space, dt; rng)
            motility = RunTumble(1.0, [30.0], 0.0)
            add_agent!(model_periodic; motility)
            adf_periodic, = run!(model_periodic, nsteps; adata)
            rng = Xoshiro(35)
            space = ContinuousSpace(extent; periodic=false)
            model_closed = StandardABM(Microbe{D}, space, dt; rng)
            add_agent!(model_closed; motility)
            adf_closed, = run!(model_closed, nsteps; adata)
            # boundary conditions don't affect vacf
            vacf_periodic = Analysis.acf(adf_periodic, :direction)
            vacf_closed = Analysis.acf(adf_closed, :direction)
            @test vacf_periodic[1] ≈ vacf_closed[1]
            @test length(vacf_periodic[1]) == nsteps+1
            # lag-0 value is unity
            @test vacf_periodic[1][1] ≈ 1
            # no value can be larger than the value at lag 0
            @test maximum(vacf_periodic[1]) == vacf_periodic[1][1]
        end
    end

    @testset "MSD" begin
        for D in 1:3
            dt = 0.1
            L = 100
            extent = fill(float(L), SVector{D})
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt)
            turn_rate = 0 # ballistic motion
            U = 30.0
            motility = RunTumble(1/turn_rate, [U], 0.0)
            add_agent!(model; motility)
            nsteps = 50
            adata = [position]
            adf, = run!(model, nsteps; adata)
            real_msd = @. (U * (0:nsteps)*dt)^2
            Analysis.unfold!(adf, model; key=:position)
            @test hasproperty(adf, :position_unfold)
            Δ² = Analysis.emsd(adf, :position_unfold)
            @test length(Δ²) == nsteps+1
            @test Δ² ≈ real_msd

            # previous should be identical in a rectangular domain
            if D > 1 # skip trivial D=1 case
                L₀ = 100; L₁ = 40
                extent = SVector{D}(i == 1 ? L₀ : L₁ for i in 1:D)
                space = ContinuousSpace(extent)
                model = StandardABM(Microbe{D}, space, dt)
                turn_rate = 0 # ballistic motion
                U = 30.0
                motility = RunTumble(1/turn_rate, [U], 0.0)
                add_agent!(model; motility)
                nsteps = 50
                adata = [position]
                adf, = run!(model, nsteps; adata)
                Analysis.unfold!(adf, model; key=:position)
                real_msd = @. (U * (0:nsteps)*dt)^2 # 30 is default speed
                Δ² = Analysis.emsd(adf, :position_unfold)
                @test Δ² ≈ real_msd
            end

            turn_rate = Inf # reversal at each step
            motility = RunReverse(1/turn_rate, [U], 1/turn_rate, [U])
            model = StandardABM(Microbe{D}, space, dt)
            add_agent!(extent./2, model; motility)
            nsteps = 3
            adata = [position]
            adf, = run!(model, nsteps; adata)
            real_msd = [0, U^2, 0, U^2] .* dt^2
            Δ² = Analysis.emsd(adf, :position) # unfolding not required here
            @test Δ² ≈ real_msd
        end
    end
end
