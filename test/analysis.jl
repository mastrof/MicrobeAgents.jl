using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

@testset "Post-processing and analysis routines" begin
    dt = 1
    L = 100
    for D in 1:3
        extent = ntuple(_ -> L, D)
        pos = ntuple(zero, D)
        U = 1
        motility = RunTumble(speed=[U])
        vel = ntuple(_ -> 1/√D, D)
        turn_rate = 0
        model = ABM(Microbe{D}, extent, dt)
        add_agent!(pos, model; vel, motility, turn_rate)
        nsteps = 10
        adata = [:pos]
        adf, mdf = run!(model, nsteps; adata)
        traj = vectorize_adf_measurement(adf, :pos)
        @test traj isa AbstractMatrix{<:Tuple}
        @test size(traj) == (nsteps+1, 1)
        real_traj = [pos .+ vel.*(U*n*dt) for n in 0:nsteps, _ in 1:1]
        Δ = norm.([traj[i] .- real_traj[i] for i in eachindex(traj)])
        for i in eachindex(Δ)
            @test Δ[i] ≈ 0 atol=1e-12
        end
    end

    @testset "Autocorrelation function" begin
        for D in 1:3
            dt = 0.1
            nsteps = 100
            v(microbe) = microbe.speed .* microbe.vel
            adata = [v]
            L = 100; extent = ntuple(_ -> L, D)
            rng = Xoshiro(35)
            model_periodic = ABM(Microbe{D}, extent, dt; rng)
            add_agent!(model_periodic)
            adf_periodic, = run!(model_periodic, nsteps; adata)
            rng = Xoshiro(35)
            model_closed = ABM(Microbe{D}, extent, dt; periodic=false, rng)
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
            turn_rate = 0 # ballistic motion
            model = ABM(Microbe{D}, extent, dt)
            add_agent!(model; turn_rate)
            nsteps = 50
            adata = [:pos]
            adf, = run!(model, nsteps; adata)
            real_msd = @. (30 * (1:nsteps-1)*dt)^2 # 30 is default speed
            Δ² = msd(adf; L) # L keyword to unfold periodic boundaries
            @test length(Δ²) == nsteps-1
            @test Δ² ≈ real_msd

            motility = RunReverse()
            turn_rate = Inf # reversal at each step
            model = ABM(Microbe{D}, extent, dt)
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