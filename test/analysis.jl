using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

@testset "Post-processing and analysis routines" begin
    dt = 1
    L = 100
    for D in 1:3
        extent = ntuple(_ -> L, D)
        space = ContinuousSpace(extent)
        pos = ntuple(zero, D)
        U = 1
        motility = RunTumble(speed=[U])
        vel = ntuple(_ -> 1/√D, D)
        turn_rate = 0
        model = StandardABM(Microbe{D}, space, dt)
        add_agent!(pos, model; vel, motility, turn_rate)
        nsteps = 10
        adata = [:pos]
        adf, mdf = run!(model, nsteps; adata)
        traj = vectorize_adf_measurement(adf, :pos)
        @test traj isa AbstractMatrix{<:SVector}
        @test size(traj) == (nsteps+1, 1)
        real_traj = [pos .+ vel.*(U*n*dt) for n in 0:nsteps, _ in 1:1]
        Δ = norm.([traj[i] .- real_traj[i] for i in eachindex(traj)])
        for i in eachindex(Δ)
            @test Δ[i] ≈ 0 atol=1e-12
        end
    end

    @testset "Turn detection and run statistics" begin
        dt = 0.1
        L = 100
        for D in 1:3
            extent = ntuple(_ -> L, D)
            space = ContinuousSpace(extent)
            model = StandardABM(Microbe{D}, space, dt)
            add_agent!(model; turn_rate=Inf) # turns every step
            nsteps = 10
            adf, _ = run!(model, nsteps; adata=[:vel])
            @test mean_turnrate(adf,dt) ≈ 1/dt
            @test mean_runduration(adf,dt) ≈ dt
            @test detect_turns(adf) == BitMatrix(ones(nsteps,1))
            @test rundurations(adf,dt) == [repeat([dt],nsteps)]

            D == 1 && continue # only test rotational diffusion for D>1
            model = StandardABM(Microbe{D}, space, dt)
            # never turns but makes small deviations due to brownian noise
            Drot = 0.1
            add_agent!(model; turn_rate=0.0, rotational_diffusivity=Drot)
            nsteps = 10
            adf, _ = run!(model, nsteps; adata=[:vel])
            # detects turns at each step if threshold is not set
            @test detect_turns(adf) == BitMatrix(ones(nsteps,1))
            # set threshold at 4σ to ignore rotational diffusion
            threshold_angle = 4 * sqrt(2*Drot*dt)
            @test detect_turns(adf;threshold_angle) == BitMatrix(zeros(nsteps,1))
            @test mean_turnrate(adf,dt;threshold_angle) == 0.0
            @test mean_runduration(adf,dt;threshold_angle) == Inf
            @test rundurations(adf,dt;threshold_angle) == [Float64[]]
        end
    end

    @testset "Autocorrelation function" begin
        for D in 1:3
            dt = 0.1
            nsteps = 100
            v(microbe) = microbe.speed .* microbe.vel
            adata = [v]
            L = 100; extent = ntuple(_ -> L, D)
            space = ContinuousSpace(extent)
            rng = Xoshiro(35)
            model_periodic = StandardABM(Microbe{D}, space, dt; rng)
            add_agent!(model_periodic)
            adf_periodic, = run!(model_periodic, nsteps; adata)
            rng = Xoshiro(35)
            space = ContinuousSpace(extent; periodic=false)
            model_closed = StandardABM(Microbe{D}, space, dt; rng)
            add_agent!(model_closed)
            adf_closed, = run!(model_closed, nsteps; adata)
            # boundary conditions don't affect vacf
            vacf_periodic = acf(adf_periodic, :v)
            vacf_closed = acf(adf_closed, :v)
            @test vacf_periodic ≈ vacf_closed
            @test length(vacf_periodic) == nsteps+1
            # lag-0 value corresponds to speed squared
            @test vacf_periodic[1] .≈ (model_periodic[1].speed)^2
            # no value can be larger than the value at lag 0
            @test maximum(vacf_periodic) == vacf_periodic[1]
        end
    end

    @testset "MSD" begin
        for D in 1:3
            dt = 0.1
            L = 100; extent = ntuple(_ -> L, D)
            space = ContinuousSpace(extent)
            turn_rate = 0 # ballistic motion
            model = StandardABM(Microbe{D}, space, dt)
            add_agent!(model; turn_rate)
            nsteps = 50
            adata = [:pos]
            adf, = run!(model, nsteps; adata)
            real_msd = @. (30 * (1:nsteps-1)*dt)^2 # 30 is default speed
            Δ² = msd(adf; L) # L keyword to unfold periodic boundaries
            @test length(Δ²) == nsteps-1
            @test Δ² ≈ real_msd

            # previous should be identical in a rectangular domain
            if D > 1 # skip trivial D=1 case
                L₀ = 100; L₁ = 40
                extent = ntuple(i -> i==1 ? L₀ : L₁, D)
                space = ContinuousSpace(extent)
                turn_rate = 0 # ballistic motion
                model = StandardABM(Microbe{D}, space, dt)
                add_agent!(model; turn_rate)
                nsteps = 50
                adata = [:pos]
                adf, = run!(model, nsteps; adata)
                real_msd = @. (30 * (1:nsteps-1)*dt)^2 # 30 is default speed
                Δ² = msd(adf; L=extent) # unfold periodic boundaries
                @test Δ² ≈ real_msd
            end

            motility = RunReverse()
            turn_rate = Inf # reversal at each step
            model = StandardABM(Microbe{D}, space, dt)
            add_agent!(extent./2, model; motility, turn_rate)
            nsteps = 5
            adata = [:pos]
            adf, = run!(model, nsteps; adata)
            real_msd = [30^2, 0, 30^2, 0] .* dt^2
            Δ² = msd(adf) # unfolding not required here
            @test Δ² ≈ real_msd
        end
    end
end
