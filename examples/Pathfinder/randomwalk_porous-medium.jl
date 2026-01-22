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
wm = walkmap(bodies, extent, R1/25, 0)
# Initialise pathfinder
pf = AStar(space; walkmap=wm)

model = StandardABM(BrownBerg{2,2}, space, dt;
    container=Vector,
    properties=Dict(:pathfinder => pf),
    agent_step! = microbe_pathfinder_step!
)
motility = RunTumble([20.0], 1.0, Isotropic2D)
for _ in 1:nbacteria
    # make sure bacteria are initialised
    # within the accessible domain space
    pos = random_position_pathfinder(model)
    add_agent!(pos, model; motility)
end
adata = [:pos]
nsteps = 1000
adf, mdf = run!(model, nsteps; adata)

traj = Analysis.adf_to_matrix(adf, :pos)
xx = range(0, extent[1], length=size(wm,1))
yy = range(0, extent[2], length=size(wm,2))
contourf(xx, yy, wm', levels=2, ratio=1, axis=false, grid=false, cbar=false, c=:bone)
plot!(first.(traj), last.(traj), leg=false)
# plot!(size=(900,900))
