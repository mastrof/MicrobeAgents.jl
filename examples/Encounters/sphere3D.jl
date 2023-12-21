using Distributed
addprocs(8; exeflags="--project")

@everywhere begin
    using Agents, MicrobeAgents
    using LinearAlgebra, Statistics
    using Random
end
using Plots

@everywhere function encounters!(model)
    for microbe in allagents(model)
        if is_encounter(microbe, model.sphere, model)
            model.encounters += 1
            reinsert!(microbe, model)
        end
        model.old_positions[microbe.id] = microbe.pos
    end
end

@everywhere function reinsert!(microbe, model)
    l = spacesize(model)[1]/2
    R = model.sphere.r + rand(abmrng(model))*(l-model.sphere.r)
    pos = l .+ rand_vel(abmrng(model), 3) .* R
    move_agent!(microbe, pos, model)
    MicrobeAgents.turn!(microbe, model)
end

function setupmodel(R, L, n; dt=0.1, rng=Xoshiro(1))
    extent = ntuple(_ -> L, 3)
    space = ContinuousSpace(extent)
    properties = Dict(
        :sphere => HyperSphere(extent./2, R),
        :encounters => 0,
    )
    model = StandardABM(Microbe{3}, space, dt; properties, rng, container=Vector)
    foreach(_ -> add_agent!(model; turn_rate=2.0), 1:n)
    model → encounters!
    abmproperties(model)[:old_positions] = map(m -> m.pos, allagents(model))
    model
end

R = Float64.([2, 5, 10, 15, 20, 25, 35, 50, 75, 100, 150, 200])
L = @. max(400,15R)
rng = map(_ -> Xoshiro(32), R)
n = 1000
model = [setupmodel(R[i], L[i], n; rng=rng[i]) for i in eachindex(R)]

@everywhere stop(model,s) = model.encounters > 2000
@everywhere when_model(model,s) = s%50 == 0
mdata = [:encounters]
mdf = (pmap(mdl -> last(run!(mdl, stop; mdata, when_model)), model))
rmprocs(workers())

t = vcat(map(m -> m.step[end] .* 0.1, mdf))
e = vcat(map(m -> m.encounters[end], mdf))
Γ = @. (e / t) * (L^3/n)
U = 30.0; τ = 0.5; D = U^2*τ/3
plot(R, Γ, scale=:log10, m=:c)
plot!(R, r -> 4π*D*r); plot!(R, r -> π*U*r^2)
