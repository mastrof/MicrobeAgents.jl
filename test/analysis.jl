using MicrobeAgents, Test
using Random
using LinearAlgebra: norm

# one ballistic microbe; returns positions and velocities every step
function ballistic_adf(D; L = 100.0, U = 1.0, dt = 1.0, nsteps = 10, vel = fill(1/√D, SVector{D}), pos = zero(SVector{D}))
    space = ContinuousSpace(fill(L, SVector{D}))
    model = StandardABM(Microbe{D}, space, dt)
    add_agent!(pos, model; vel, motility = RunTumble([U], Inf, Isotropic(D)))
    adf, = run!(model, nsteps; adata = [position, velocity])
    return adf, model
end

@testset "Analysis" begin
    @testset "Dataframe conversion" begin
        for D in 1:3
            nsteps = 10
            vel = fill(1/√D, SVector{D})
            adf, = ballistic_adf(D; nsteps, vel)
            traj = Analysis.adf_to_matrix(adf, :position)
            @test traj isa AbstractMatrix{<:SVector}
            @test size(traj) == (nsteps+1, 1)
            @test traj[:, 1] ≈ [vel .* n for n in 0:nsteps]
            @test Analysis.adf_to_vectors(adf, :position) == [traj[:, 1]]
        end
    end

    @testset "Unfolding periodic trajectories" begin
        for D in 1:3
            # crosses the periodic boundary several times
            L, U, nsteps = 10.0, 3.0, 20
            vel = fill(1/√D, SVector{D})
            adf, model = ballistic_adf(D; L, U, nsteps, vel)
            Analysis.unfold!(adf, model; key = :position)
            @test hasproperty(adf, :position_unfold)
            @test adf.position_unfold ≈ [vel .* (U * n) for n in 0:nsteps]
        end
    end

    @testset "Turn detection" begin
        space = ContinuousSpace((100.0, 100.0))
        nsteps = 10
        model = StandardABM(Microbe{2}, space, 1.0; rng = Xoshiro(1))
        add_agent!(model; motility = RunReverse([1.0], 0.0, [1.0], 0.0)) # reverses every step
        add_agent!(model; motility = RunTumble([1.0], Inf, Isotropic(2))) # never turns
        add_agent!(model; motility = RunTumble([1.0], Inf, Isotropic(2)), rotational_diffusivity = 0.01)
        adf, = run!(model, nsteps; adata = [velocity])
        Analysis.detect_turns!(adf)
        turns(i, key = :has_turned) = adf[adf.id .== i, key]
        @test turns(1) == [false; trues(nsteps)]
        @test turns(2) == falses(nsteps + 1)
        @test all(turns(3)[2:end]) # any rotation counts without a threshold
        @test Analysis.run_durations(adf)[1:2] == [ones(Int, nsteps), Int[]]
        # a threshold well above the rotational noise ignores it
        Analysis.detect_turns!(adf; threshold_angle = 4 * sqrt(2 * 0.01 * 1.0), new_key = :big_turn)
        @test !any(turns(3, :big_turn))
        @test Analysis.run_durations(adf; turns_key = :big_turn)[1] == ones(Int, nsteps)
    end

    @testset "Drift velocity" begin
        U = 2.0
        adf, = ballistic_adf(2; U, vel = SVector(1.0, 0.0), pos = SVector(10.0, 50.0))
        Analysis.driftvelocity_direction!(adf, (1.0, 0.0))
        @test all(adf.drift_direction .≈ U)
        Analysis.driftvelocity_direction!(adf, SVector(0.0, 3.0))
        @test all(abs.(adf.drift_direction) .< 1e-12)
        Analysis.driftvelocity_direction!(adf, SVector(1.0, 1.0); normalize = true)
        @test all(adf.drift_direction .≈ 1/√2)
        # towards a point straight ahead / behind
        Analysis.driftvelocity_point!(adf, (90.0, 50.0))
        @test all(adf.drift_point .≈ U)
        Analysis.driftvelocity_point!(adf, SVector(0.0, 50.0); normalize = true)
        @test all(adf.drift_point .≈ -1)
    end
end
