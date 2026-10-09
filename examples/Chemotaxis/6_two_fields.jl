# # Two opposing chemical fields

#=
Microbes can respond to several chemicals at once. Each chemical is a model
property holding a `ChemicalField`; `:chemicalfield` is the default one.
A behavior senses the field named by its `field` keyword, so a microbe that
is attracted to one substance and repelled by another carries one behavior
per substance, each with its own parameters.

Here the attractant increases with `x` while the repellent decreases with `x`,
so that both gradients push towards positive `x`.
We compare a population sensing only the attractant with one sensing both.
=#
using MicrobeAgents
using Plots
using Random

L = 6000.0 # μm
space = ContinuousSpace((L, L); periodic = false)
attractant = ChemicalField{2}(;
    concentration_field = (m, model) -> position(m)[1] / L * 10,
    concentration_gradient = (m, model) -> SVector(10 / L, 0.0),
)
repellent = ChemicalField{2}(;
    concentration_field = (m, model) -> (L - position(m)[1]) / L * 10,
    concentration_gradient = (m, model) -> SVector(-10 / L, 0.0),
)
dt = 0.1 # s
model = StandardABM(Microbe{2}, space, dt;
    properties = Dict(:chemicalfield => attractant, :repellent => repellent),
    rng = Xoshiro(42),
)

#=
The attractant-only population carries one `BrownBerg`; the other carries a second
one tied to the `:repellent` field, with a negative `gain`.
Biases multiply across behaviors, so the two signals combine.
=#
motility = RunTumble([30.0], 1.0, Isotropic(2))
attract = BrownBerg(gain = 660, receptor_binding_constant = 5)
repel = BrownBerg(field = :repellent, gain = -660, receptor_binding_constant = 5)
center = SVector(L / 2, L / 2)
for _ in 1:100
    add_agent!(center, model; motility, behaviors = (attract,))        # attractant only
    add_agent!(center, model; motility, behaviors = (attract, repel))  # attractant + repellent
end

#=
Agents are added alternately, so odd ids belong to the attractant-only population
and even ids to the other one.
We record the mean `x` of each population along the run.
=#
meanx(ids) = sum(position(model[i])[1] for i in ids) / length(ids)
attr_ids = 1:2:200
both_ids = 2:2:200

nsteps = 1000 # 100 s, before the front reaches the wall
x_attr = [meanx(attr_ids)]
x_both = [meanx(both_ids)]
for _ in 1:nsteps
    run!(model, 1)
    push!(x_attr, meanx(attr_ids))
    push!(x_both, meanx(both_ids))
end

t = (0:nsteps) .* dt
plot(t, x_attr; lw = 2, label = "attractant")
plot!(t, x_both; lw = 2, label = "attractant + repellent")
plot!(xlabel = "time (s)", ylabel = "mean x (μm)", legend = :topleft)

#=
Both populations climb the attractant gradient, and the second one is also
pushed away from the repellent, so it drifts faster.
=#
