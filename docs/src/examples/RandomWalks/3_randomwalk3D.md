```@meta
EditURL = "../../../../examples/RandomWalks/3_randomwalk3D.jl"
```

# 3D Random walk

Without any significant difference, we can also simulate three-dimensional random walks.

````@example 3_randomwalk3D
using MicrobeAgents
using Distributions
using Plots

L = 500
space = ContinuousSpace((L,L,L))
dt = 0.1
model = StandardABM(Microbe{3}, space, dt)

add_agent!(model; motility=RunReverse(1.0, [55], 1.0, [55]), rotational_diffusivity=0.2)
add_agent!(model; motility=RunTumble(2.0, Normal(30,6), Isotropic(3)))
add_agent!(model; motility=RunReverseFlick(1.0, [30], 1.0, [6]), rotational_diffusivity=0.1)

nsteps = 600
adata = [position]
adf, _ = run!(model, nsteps; adata)

Analysis.unfold!(adf, model)
traj = Analysis.adf_to_matrix(adf, :position_unfold)
x = first.(traj)
y = getindex.(traj, 2)
z = last.(traj)
t = axes(x,1) .* dt

plot(x, y, z, xlab="x", ylab="y", zlab="z", lw=2, ratio=1,
    lab=["RunReverse" "RunTumble" "RunReverseFlick"]
)
````

