```@meta
EditURL = "../../../../examples/Chemotaxis/4_response_functions.jl"
```

# Comparison of chemotactic response functions

Here we will compare the chemotactic response function of the `Celani`
and `BrownBerg` model to an impulse stimulus of chemoattractant.

While `Celani` only needs the `concentration_field` to determine the
chemotactic response, `BrownBerg` also needs the the time derivative (`concentration_ramp`)
to be defined explicitly (also the `concentration_gradient` but it's not relevant in this specific study).

````@example 4_response_functions
using MicrobeAgents
using Plots

θ(a,b) = a>b ? 1.0 : 0.0 # Heaviside theta function
function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    t₁ = model.t₁
    t₂ = model.t₂
    dt = model.timestep
    t = abmtime(model) * dt
    concentration_field(t, C₀, C₁, t₁, t₂)
end
concentration_field(t,C₀,C₁,t₁,t₂) = C₀+C₁*θ(t,t₁)*(1-θ(t,t₂))

δ(t,dt) = 0 <= t <= dt ? 1.0/dt : 0.0 # discrete approximation to Dirac delta
function concentration_ramp(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    t₁ = model.t₁
    t₂ = model.t₂
    dt = model.timestep
    t = abmtime(model) * dt
    concentration_ramp(t, C₀, C₁, t₁, t₂, dt)
end
function concentration_ramp(t, C₀, C₁, t₁, t₂, dt)
    C₁*(δ(t-t₁, dt) - δ(t-t₂, dt))
end

space = ContinuousSpace(ntuple(_ -> 500.0, 3)) # μm
C₀ = 1.0 # μM
C₁ = 2.0-C₀ # μM
T = 50.0 # s
dt = 0.1 # s
t₁ = 10.0 # s
t₂ = 30.0 # s
properties = Dict(
    :chemoattractant => GenericChemoattractant{3,Float64}(;
        concentration_field,
        concentration_ramp
    ),
    :C₀ => C₀,
    :C₁ => C₁,
    :t₁ => t₁,
    :t₂ => t₂,
)

model = StandardABM(Union{BrownBerg{3},Celani{3}}, space, dt; properties)

add_agent!(BrownBerg{3}, model; motility=RunTumble([0], 1.0, Isotropic(3)),
    memory=1,
)
add_agent!(Celani{3}, model; motility=RunTumble([0], 1.0, Isotropic(3)), gain=4)
add_agent!(Celani{3}, model; motility=RunTumble([0], 1.0, Isotropic(3)), gain=4,
    chemotactic_precision=50.0
)

nsteps = round(Int, T/dt)
adata = [bias]
adf, = run!(model, nsteps; adata)

S = Analysis.adf_to_matrix(adf, :bias)

_pink = palette(:default)[4]
plot()
x = (0:dt:T) .- t₁
plot!(
    x, S,
    lw=1.5, lab=["BrownBerg" "Celani" "Celani + Noise"]
)
plot!(ylims=(-0.1,2.1), ylab="Response", xlab="time (s)")
plot!(twinx(),
    x, t -> concentration_field(t.+t₁,C₀,C₁,t₁,t₂),
    ls=:dash, lw=1.5, lc=_pink, lab=false,
    tickfontcolor=_pink,
    ylab="C (μM)", guidefontcolor=_pink
)
````

