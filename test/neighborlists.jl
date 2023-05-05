using Test
using MicrobeAgents
using Agents
using StaticArrays

@testset "Neighbor Lists" begin
    @testset "StandardABM" begin
        for D in 2:3
            L = 100
            extent = ntuple(_->L, D)
            dt = 1
            model = StandardABM(Microbe{D}, extent, dt)
            n = 10
            foreach(_ -> add_agent!(model), 1:n)
            listkey = :neighbors
            cutoff = 10.0
            neighborlist!(model, cutoff, listkey)
            # check that the key is correctly added to the model
            @test haskey(abmproperties(model), listkey)
            # check that positions are correctly stored in the neighbor list
            microbe_positions = [SVector(model[i].pos) for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # add neighbor list updater to model
            model → (model) -> update_neighborlist!(model, listkey)
            run!(model)
            # after the step, the neighbor list positions should be updated to the new positions
            microbe_positions = [SVector(model[i].pos) for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # if a microbe is removed its position is not updated anymore and ids are retained
            j = 6
            remove_agent!(j, model)
            pos_removed = microbe_positions[j]
            run!(model, 20)
            microbe_positions = [i==j ? pos_removed : SVector(model[i].pos) for i in 1:n]
            @test microbe_positions == model.neighbors.xpositions
        end
    end

    @testset "UnremovableABM" begin
        for D in 2:3
            L = 100
            extent = ntuple(_->L, D)
            dt = 1
            model = UnremovableABM(Microbe{D}, extent, dt)
            n = 10
            foreach(_ -> add_agent!(model), 1:n)
            listkey = :neighbors
            cutoff = 10.0
            neighborlist!(model, cutoff, listkey)
            # check that the key is correctly added to the model
            @test haskey(abmproperties(model), listkey)
            # check that positions are correctly stored in the neighbor list
            microbe_positions = [SVector(model[i].pos) for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # add neighbor list updater to model
            model → (model) -> update_neighborlist!(model, listkey)
            run!(model)
            # after the step, the neighbor list positions should be updated to the new positions
            microbe_positions = [SVector(model[i].pos) for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
        end
    end
end
