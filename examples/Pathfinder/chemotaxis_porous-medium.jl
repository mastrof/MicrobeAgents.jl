using MicrobeAgents
using DelimitedFiles
using BubbleBath
using LinearAlgebra
using Plots

# Concentration field from left to right
function concentration_field(pos, model)
    x,y = pos
    Lx,Ly = spacesize(model)
    C₁, C₂ = model.C₁, model.C₂
    concentration_field(x,y,Lx,Ly,C₁,C₂)
end
concentration_field(x,y,Lx,Ly,C₁,C₂) = C₁ + x/Lx*(C₂-C₁)
function concentration_gradient(pos, model)
    x,y = pos
    Lx,Ly = spacesize(model)
    C₁, C₂ = model.C₁, model.C₂
    ((C₂-C₁)/Lx, 0.0)
end

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

C₁, C₂ = 0.0, 10.0
properties = Dict(
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient,
    :C₁ => C₁,
    :C₂ => C₂,
)

model = StandardABM(BrownBerg{2}, space, dt; properties, container=Vector)
pathfinder!(model, wm)
foreach(_ -> add_agent!((0.0, rand()*extent[2]), model), 1:nbacteria)
adata = [:pos]
nsteps = 5000
adf, mdf = run!(model,
    microbe_pathfinder_step!, model.update!, nsteps;
    adata
)

traj = vectorize_adf_measurement(adf, :pos)
xx = range(0, extent[1], length=size(wm,1))
yy = range(0, extent[2], length=size(wm,2))
cxy = concentration_field.(xx',yy,extent...,C₁,C₂)
contourf(xx, yy, cxy,
    cbar=false, levels=30, ratio=1, axis=false, grid=false,
    lw=0, c=:viridis, bgcolor=:black, size=(800,400),
)
zz = Float64.(wm')
zz[zz.==1] .= NaN
contourf!(xx, yy, zz, levels=1, ratio=1, axis=false, grid=false, cbar=false, c=:bone)
plot!(first.(traj), last.(traj), leg=false, palette=:Reds)

#=
x = first.(traj); y = last.(traj)
for t in range(2,nsteps+1,step=10)
    contourf(xx, yy, cxy,
        cbar=false, levels=30, ratio=1, axis=false, grid=false,
        lw=0, c=:viridis, bgcolor=:black, size=(800,400),
    )
    contourf!(xx,yy,zz,levels=1,cbar=false,c=:bone)
    xt = @view x[1:t,:]
    yt = @view y[1:t,:]
    plot!(xt, yt, leg=false, palette=:Reds)
    it = lpad(t, 4, '0')
    savefig("chemo_anim/frame_$(it).png")
end
=#
