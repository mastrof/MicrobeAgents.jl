# Chemotaxis

## Linear concentration ramp
As a first step into the world of chemotaxis, we will reproduce
an in-silico version of a classical laboratory assay: a linear concentration
ramp in a rectangular channel. We will use the `BrownBerg` model.

To define concentration profiles, MicrobeAgents.jl expects three functions:
`concentration_field`, `concentration_gradient` and `concentration_time_derivative`,
all of which must have a method with signature `(pos,model)`, i.e. they take
as input a microbe position and the ABM object.
In this example the field is static, so we won't define a time derivative.

```
using MicrobeAgents

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
```
The quantities `C₀` and `C₁` represent the concentration values at left
(`x=0`) and right (`x=L`) edges of the channel.
Notice that `concentration_gradient` returns a `Tuple` with the same
size as `pos`, not just a `Number`.
This is required since the gradient is a vector quantity.

We can set up the ABM as usual, but we will need to supply `:concentration_field`
and `:concentration_gradient` to the `properties` container, as well as
`:C₀` and `:C₁`. Then we can run the simulation and make a nice plot.
```
Lx, Ly = 3000, 1500
extent = (Lx, Ly)
periodic = false
space = ContinuousSpace(extent; periodic)
Δt = 0.1 # timestep (s)
T = 120 # simulation time (s)
nsteps = round(Int, T/Δt)
n = 100

C₀, C₁ = 0.0, 20.0 # μM
properties = Dict(
    :C₀ => C₀,
    :C₁ => C₁,
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient
)
model = StandardABM(BrownBerg{2}, space, Δt; periodic, properties, container=Vector)
for i in 1:n
    add_agent!(model)
end

adata = [:pos]
adf, mdf = run!(model, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
y = last.(traj)

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
```
![Drift of chemotactic bacteria in a linear concentration ramp](linear-ramp.svg)
