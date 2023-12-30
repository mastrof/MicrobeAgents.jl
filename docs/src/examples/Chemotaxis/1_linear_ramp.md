```@meta
EditURL = "../../../../examples/Chemotaxis/1_linear_ramp.jl"
```

# Linear concentration ramp

In this example we will setup an in-silico version of a typical laboratory assay,
with chemotactic bacteria moving in a linear concentration ramp, i.e. a concentration
field of the form
```math
C(x) = C_0 + (C_1 - C_0)\dfrac{x}{L_x}, \qquad x \in [0,L_x].
```

We will create a closed two-dimensional domain, mimicking a thin microfluidic chamber,
with a length of 3 mm along the `x` direction, and 1.5 mm along the `y` direction.
We will then setup the concentration field along the `x` direction and observe
chemotactic microbes as they drift towards the high-concentration region of the chamber.

The first thing we have to do is define two functions for the `concentration_field` and the
`concentration_gradient`. They must take as arguments the position of a single microbe,
and the model (from which we can access other properties of the system).
Of course, for our convenience we can dispatch these functions on whatever arguments we want,
as long as they have a method whose signature matches the MicrobeAgents API.

Importantly, the `concentration_field` must return a scalar, non-negative value.
Since the gradient is a vector quantity, the `concentration_gradient` should instead return
an iterable with length equal to the system dimensionality; a `SVector` is the recommended
choice, but `NTuple`s, `Vector`s, etc.. work just fine.

All the parameters that we need to evaluate the concentration field and gradient, in our case
the two concentration values `C₀` and `C₁` and the chamber length `Lx`, should be extracted
from the `model`.

````@example 1_linear_ramp
using MicrobeAgents
using Plots

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
@inline concentration_gradient(pos,Lx,C₀,C₁) = SVector{length(pos)}(i==1 ? (C₁-C₀)/Lx : 0.0 for i in eachindex(pos))
````

Now as usual we define the simulation domain and the integration timestep, but we also define
a `properties` dictionary, which we pass as a keyword argument to `StandardABM`.
This dictionary will contain all the information regarding our concentration field.

Note that the `:C₀` and `:C₁` keys have been defined by us; we could have chosen different
names for them.
The `concentration_field` and `concentration_gradient` functions instead **must** be assigned
to the `:concentration_field` and `:concentration_gradient` keys respectively; this is required
by MicrobeAgents and assigning these functions to any other key will not produce the desired results.

To observe chemotaxis, we must use a microbe type for which chemotactic behavior is implemented.
If we used the base `Microbe`, no matter what field we define, we would only observe a random walk
since no chemotactic behavior is implemented.
The most classic model of chemotaxis is implemented in the `BrownBerg` type;
we will not modify its parameters here and just stick to the default values.

````@example 1_linear_ramp
Lx, Ly = 3000, 1500 # domain size (μm)
periodic = false
space = ContinuousSpace((Lx,Ly); periodic)
Δt = 0.1 # timestep (s)

# model setup
C₀, C₁ = 0.0, 20.0 # μM
properties = Dict(
    :C₀ => C₀,
    :C₁ => C₁,
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient
)
model = StandardABM(BrownBerg{2}, space, Δt; properties)
n = 100 # number of microbes
for i in 1:n
    add_agent!(model) # add at random positions
end
model
````

Now that the model is created, we just run it as usual, collecting the position
of the microbes at each timestep.
The visualization is slightly more involved since we want to plot
the microbe trajectories on top of the concentration field shown as a heatmap,
but there is really no difference from what we have seen in the random walk examples.

In the figure, we will see that all the microbes drift towards the right,
where the concentration of the attractant is higher.

````@example 1_linear_ramp
T = 120 # simulation time (s)
nsteps = round(Int, T/Δt)
adata = [position]
adf, mdf = run!(model, nsteps; adata)

traj = Analysis.adf_to_matrix(adf, :position)
x = first.(traj)
y = last.(traj)

ts = unique(adf.step) .* Δt
lw = eachindex(ts) ./ length(ts) .* 3
xmesh = range(0,Lx,length=100)
ymesh = range(0,Ly,length=100)
xn = @view x[:,1:10]
yn = @view y[:,1:10]
c = concentration_field.(Iterators.product(xmesh,ymesh),Lx,C₀,C₁)
heatmap(xmesh, ymesh, c', cbar=false, c=:bone,
    ratio=1, axis=false, grid=false, xlims=(0,Lx), ylims=(0,Ly)
)
plot!(xn, yn, lab=false, lw=lw, lc=(1:n)')
scatter!(xn[end,:], yn[end,:], lab=false, m=:c, mc=1:n, msw=0.5, ms=8)
````

