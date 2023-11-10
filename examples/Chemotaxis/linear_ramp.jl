using MicrobeAgents
using Plots

## define concentration field
@inline function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(spacesize(model))
    concentration_field(pos,Lx,C₀,C₁)
end
@inline concentration_field(pos,Lx,C₀,C₁) = C₀ + (C₁-C₀)*pos[1]/Lx

@inline function concentration_gradient(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(spacesize(model))
    concentration_gradient(pos,Lx,C₀,C₁)
end
@inline concentration_gradient(pos,Lx,C₀,C₁) = ntuple(i->i==1 ? (C₁-C₀)/Lx : 0.0, length(pos))

## simulation parameters
Lx, Ly = 3000, 1500 # domain size (μm)
periodic = false
space = ContinuousSpace((Lx,Ly); periodic)
Δt = 0.1 # timestep (s)
T = 120 # simulation time (s)
nsteps = round(Int, T/Δt)
n = 100

## model setup
C₀, C₁ = 0.0, 20.0 # μM
properties = Dict(
    :C₀ => C₀,
    :C₁ => C₁,
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient
)
model = StandardABM(BrownBerg{2}, space, Δt; properties)
for i in 1:n
    add_agent!(model)
end

## run
adata = [:pos, :vel]
adf, mdf = run!(model, nsteps; adata)

## postprocessing
traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
y = last.(traj)

## plotting
ts = unique(adf.step) .* Δt
lw = eachindex(ts) ./ length(ts) .* 3
xmesh = range(0,Lx,length=100)
ymesh = range(0,Ly,length=100)
xn = @view x[:,1:10]
yn = @view y[:,1:10]
c = concentration_field.(Iterators.product(xmesh,ymesh),Lx,C₀,C₁)
heatmap(xmesh, ymesh, c', cbar=false, ratio=1, axis=false, c=:bone)
plot!(xn, yn, lab=false, lw=lw, lc=(1:n)')
scatter!(xn[end,:], yn[end,:], lab=false, m=:c, mc=1:n, msw=0.5, ms=8)
plot!(size=(600,300), margin=-30Plots.px)
