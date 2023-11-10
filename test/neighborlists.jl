using Test
using MicrobeAgents
using Agents

@testset "Neighbor Lists" begin
    @testset "StandardABM" begin
        for D in 1:3
            L = 100
            extent = fill(float(L), SVector{D})
            space = ContinuousSpace(extent)
            dt = 1
            model = StandardABM(Microbe{D}, space, dt)
            n = 10
            foreach(_ -> add_agent!(model), 1:n)
            listkey = :neighbors
            cutoff = 10.0
            if D == 1
                @test_throws ArgumentError neighborlist!(model, cutoff, listkey)
                continue
            end
            neighborlist!(model, cutoff, listkey)
            # check that the key is correctly added to the model
            @test haskey(abmproperties(model), listkey)
            # check that positions are correctly stored in the neighbor list
            microbe_positions = [model[i].pos for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # add neighbor list updater to model
            model → (model) -> update_neighborlist!(model, listkey)
            run!(model)
            # after the step, the neighbor list positions should be updated to the new positions
            microbe_positions = [model[i].pos for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # if a microbe is removed its position is not updated anymore and ids are retained
            j = 6
            remove_agent!(j, model)
            pos_removed = microbe_positions[j]
            run!(model, 20)
            microbe_positions = [i==j ? pos_removed : model[i].pos for i in 1:n]
            @test microbe_positions == model.neighbors.xpositions
        end
    end

    @testset "UnremovableABM" begin
        for D in 1:3
            L = 100
            extent = fill(float(L), SVector{D})
            space = ContinuousSpace(extent)
            dt = 1
            model = StandardABM(Microbe{D}, space, dt; container=Vector)
            n = 10
            foreach(_ -> add_agent!(model), 1:n)
            listkey = :neighbors
            cutoff = 10.0
            if D == 1
                @test_throws ArgumentError neighborlist!(model, cutoff, listkey)
                continue
            end
            neighborlist!(model, cutoff, listkey)
            # check that the key is correctly added to the model
            @test haskey(abmproperties(model), listkey)
            # check that positions are correctly stored in the neighbor list
            microbe_positions = [model[i].pos for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
            # add neighbor list updater to model
            model → (model) -> update_neighborlist!(model, listkey)
            run!(model)
            # after the step, the neighbor list positions should be updated to the new positions
            microbe_positions = [model[i].pos for i in sort(collect(allids(model)))]
            @test microbe_positions == model.neighbors.xpositions
        end
    end
end
