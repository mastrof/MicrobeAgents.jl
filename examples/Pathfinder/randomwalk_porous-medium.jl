using MicrobeAgents
using DelimitedFiles
using BubbleBath
using LinearAlgebra
using Plots

# Physical parameters
dt = 0.1 # s
extent = (1000.0, 500.0) # μm
space = ContinuousSpace(extent; periodic=false)
periodic = false
nbacteria = 10

# Initialise obstacles from file
obstacle_data = readdlm("phi065_rmin5_Lx1000_Ly500.dat")
bodyradii = obstacle_data[:,1] # μm
min_radius = minimum(bodyradii)
bodypositions = [Tuple(obstacle_data[i,2:3]) for i in axes(obstacle_data,1)]
bodies = [
    Sphere(pos,r) for (r,pos) in zip(bodyradii, bodypositions)
]
wm = walkmap(bodies, extent, min_radius/25, 0)

model = StandardABM(BrownBerg{2}, space, dt; container=Vector)
pathfinder!(model, wm)
foreach(_ -> add_agent!(model), 1:nbacteria)
adata = [:pos]
nsteps = 2000
adf, mdf = run!(model, microbe_pathfinder_step!, model.update!, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
xx = range(0, extent[1], length=size(wm,1))
yy = range(0, extent[2], length=size(wm,2))
contourf(xx, yy, wm', levels=2, ratio=1, axis=false, grid=false, cbar=false, c=:bone)
plot!(first.(traj), last.(traj), leg=false)
