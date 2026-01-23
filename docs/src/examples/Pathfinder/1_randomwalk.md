```@meta
EditURL = "../../../../examples/Pathfinder/1_randomwalk.jl"
```

# Random walk in porous medium

MicrobeAgents can interface with the Agents.Pathfinding module
which allows you to perform simulations in spatially
structured environments.

We generate a simple representation of a porous medium
using the BubbleBath package, which produces a random packing
of spherical obstacles with user-specified packing fraction
and size distribution.
The set of positions and radii of these spheres are converted
into a "walkmap", a boolean array where `true` represents areas
that are accessible to the bacteria, and `false` identifies
inaccessible areas.

When the pathfinder is included in the model properties,
and a pathfinder-aware stepping function is used
(MicrobeAgents provides `microbe_pathfinder_step!` that is
fully equivalent to `microbe_step!` with added pathfinding),
the bacteria will perform a random walk as usual but
will halt when hitting the spherical obstacles.

**Important**: the pathfinder does not support continuous
geometries, only discrete representations of them.
If smooth geometric features are crucial to your application,
this approach is not suitable.

````@example 1_randomwalk
using MicrobeAgents
using DelimitedFiles
using Distributions
using BubbleBath
using LinearAlgebra
using Plots

# Physical parameters
dt = 0.1 # s
extent = (500.0, 500.0) # μm
space = ContinuousSpace(extent; periodic=false)
periodic = false
nbacteria = 10

# Generate random spherical obstacles using BubbleBath
R1, R2 = 10.0, 20.0 # min and max obstacle radii
φ = 0.3 # packing fraction
bodies = bubblebath(Uniform(R1, R2), φ, extent)
resolution = R1/25
wm = walkmap(bodies, extent, resolution, 0)

# Initialise pathfinder and add it to model properties
pf = AStar(space; walkmap=wm)
model = StandardABM(BrownBerg{2,2}, space, dt; container=Vector,
    # the key *must* be :pathfinder
    properties=Dict(:pathfinder => pf),
    # use the pathfinder-aware stepping function
    agent_step! = microbe_pathfinder_step!
)
# populate model
motility = RunTumble([20.0], 1.0, Isotropic2D)
for _ in 1:nbacteria
    # random_position_pathfinder samples random positions
    # in the model domain only within the accessible area
    pos = random_position_pathfinder(model)
    add_agent!(pos, model; motility)
end
adata = [position]
nsteps = 1000
adf, mdf = run!(model, nsteps; adata)

traj = Analysis.adf_to_matrix(adf, :position)
xx = range(0, extent[1], length=size(wm,1))
yy = range(0, extent[2], length=size(wm,2))
contourf(xx, yy, wm', levels=2, ratio=1, axis=false, grid=false, cbar=false, c=:bone)
plot!(first.(traj), last.(traj), leg=false)
````

