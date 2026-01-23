# # Linear concentration ramp with obstacle

#=
Since the pathfinder embeds with all the other funcitonalities
of the library, we can now replicate the chemotaxis example
with the linear ramp, with the addition of some obstacles:
We will add a vertical wall that prevents horizontal drift
along the gradient direction, with only a little hole for
bacteria to pass through.
=#

using MicrobeAgents
using Plots

@inline function concentration_field(microbe, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(spacesize(model))
    pos = position(microbe)
    concentration_field(pos,Lx,C₀,C₁)
end
@inline concentration_field(pos,Lx,C₀,C₁) = C₀ + (C₁-C₀)*pos[1]/Lx

@inline function concentration_gradient(microbe, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(spacesize(model))
    pos = position(microbe)
    concentration_gradient(pos,Lx,C₀,C₁)
end
@inline concentration_gradient(pos,Lx,C₀,C₁) = SVector{length(pos)}(i==1 ? (C₁-C₀)/Lx : 0.0 for i in eachindex(pos))

Lx, Ly = 3000, 1500 # domain size (μm)
periodic = false
space = ContinuousSpace((Lx,Ly); periodic)
Δt = 0.05 # timestep (s)

## model setup
C₀, C₁ = 0.0, 20.0 # μM
wm = BitMatrix(ones(Lx, Ly))
wm[1400:1600, 1:650] .= false
wm[1400:1600, 850:1500] .= false
pathfinder = AStar(space; walkmap=wm)
properties = Dict(
    :C₀ => C₀,
    :C₁ => C₁,
    :chemoattractant => GenericChemoattractant{2}(;
        concentration_field, concentration_gradient
    ),
    :pathfinder => pathfinder
)
model = StandardABM(BrownBerg{2,2}, space, Δt;
    properties,
    agent_step! = microbe_pathfinder_step!
)
n = 20 # number of microbes
for i in 1:n
    ## start all bacteria in left half of the domain
    x = rand(abmrng(model)) * 700
    y = rand(abmrng(model)) * Ly
    pos = SVector{2,Float64}(x,y)
    motility = RunTumble([30.0], 0.67, Isotropic(2))
    rotational_diffusivity = 0.1
    add_agent!(pos, model; motility, rotational_diffusivity)
end
model

T = 300 # simulation time (s)
nsteps = round(Int, T/Δt)
adata = [position]
adf, mdf = run!(model, nsteps; adata, when=5)

traj = Analysis.adf_to_matrix(adf, :position)
x = first.(traj)
y = last.(traj)

ts = unique(adf.time) .* Δt
lw = eachindex(ts) ./ length(ts) .* 3
xmesh = range(0,Lx,length=100)
ymesh = range(0,Ly,length=100)
c = concentration_field.(Iterators.product(xmesh,ymesh),Lx,C₀,C₁)
heatmap(xmesh, ymesh, c', cbar=false, c=:Greens,
    ratio=1, axis=false, grid=false, xlims=(0,Lx), ylims=(0,Ly)
)
onetonan(x) = isone(x) ? NaN : float(x)
contourf!(onetonan.(wm)', cbar=false, c=:black)
plot!(x, y, lab=false, lw=lw, lc=(1:n)', palette=palette(:acton,size(x,2)))
scatter!(x[end,:], y[end,:], lab=false, m=:c, mc=1:n, msw=0.5, ms=8)
