```@meta
EditURL = "../../../../examples/Chemotaxis/2_celani_gauss2D.jl"
```

# Noisy chemotaxis towards Gaussian source

In this example we set up a static Gaussian source and observe the chemotactic behavior
of the `Celani` model, in the presence of sensing noise (via the `chemotactic_precision`).
Playing with the `chemotactic_precision`, it can be seen that the clustering of bacteria
at the source becomes stronger with decreasing noise (decreasing chemotactic precision).

````@example 2_celani_gauss2D
using MicrobeAgents
using Plots

function concentration_field(pos, model)
    C = model.C
    σ = model.σ
    p₀ = model.p₀
    concentration_field(pos, p₀, C, σ)
end
concentration_field(pos, p₀, C, σ) = C * exp(-sum(abs2.(pos.-p₀))/(2*σ^2))

timestep = 0.1 # s
extent = ntuple(_ -> 1000.0, 2) # μm
space = ContinuousSpace(extent; periodic=false)
p₀ = extent./2 # μm
C = 20.0 # μM
σ = 100.0 # μm
properties = Dict(
    :chemoattractant => GenericChemoattractant{2,Float64}(; concentration_field),
    :C => C,
    :σ => σ,
    :p₀ => p₀,
)

model = StandardABM(Celani{2,2}, space, timestep; properties)
motility = RunTumble([30.0], 0.67, Isotropic(2); tumble_duration=0.1)

for _ in 1:300
    add_agent!(model; motility,
        chemotactic_precision=6.0, rotational_diffusivity=0.1
    )
end

nsteps = 600
adata = [position, bias]
adf, = run!(model, nsteps; adata)

traj = Analysis.adf_to_matrix(adf, :position)
xmesh = range(0, first(spacesize(model)); length=80)
ymesh = range(0, last(spacesize(model)); length=80)
c = [concentration_field(p, p₀, C, σ) for p in Iterators.product(xmesh, ymesh)]
heatmap(xmesh, ymesh, c', cbar=false, ratio=1, c=:bone, axis=false)
x = getindex.(traj,1)[end-400:5:end, :]
y = getindex.(traj,2)[end-400:5:end, :]
a = axes(x,1) ./ size(x,1)
plot!(x, y,
    lab=false, lims=(0,1000), lw=a.^2, alpha=a,
)
````

