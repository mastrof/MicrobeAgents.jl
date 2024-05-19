```@meta
EditURL = "../../../../examples/Chemotaxis/3_xie_response-function.jl"
```

# Response function (Xie)

In this example we will probe the response function implemented
in the `Xie` model of chemotaxis.
The impulse response function is the "output" of the bacterial chemotaxis
pathway when presented with an input signal`.

To do this, we will emulate another classical laboratory assay, where the
bacterium is tethered to a wall, and it is exposed to a temporal change
in the concentration of a chemoattractant.
The response to the stimulus can be measured by observing modulations
in the instantaneous tumbling rate.
For each of the implemented microbe types, MicrobeAgents provides a
`bias` function which returns the instantaneous bias
in the tumbling rate, evaluated from the internal state of the microbe.
Monitoring the time evolution of the tumble bias under teporal stimuli
then allows us to access the response function of the microbe.

In the `Xie` model, chemotaxis is implemented by direct samplings of the
`concentration_field`, thus we don't need to explicitly define neither
a `concentration_gradient` nor a `concentration_time_derivative`.
We will represent our temporal stimuli in the form of square waves which
instantaneously switch from a baseline value `C₀` to a peak value `C₁+C₀`
homogeneously over space.
The excitation will occur at a time `t₁` and go back to baseline levels
at a time `t₂`.

````@example 3_xie_response-function
using MicrobeAgents
using Plots

θ(a,b) = a>b ? 1.0 : 0.0 # heaviside theta function
function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    t₁ = model.t₁
    t₂ = model.t₂
    dt = model.timestep
    t = abmtime(model) * dt
    # notice the time dependence!
    concentration_field(t, C₀, C₁, t₁, t₂)
end
concentration_field(t,C₀,C₁,t₁,t₂) = C₀+C₁*θ(t,t₁)*(1-θ(t,t₂))

space = ContinuousSpace(ntuple(_ -> 500.0, 3)) # μm
C₀ = 0.01 # μM
C₁ = 5.0-C₀ # μM
T = 50.0 # s
t₁ = 15.0 # s
t₂ = 35.0 # s
properties = Dict(
    :chemoattractant => GenericChemoattractant{3,Float64}(; concentration_field),
    :C₀ => C₀,
    :C₁ => C₁,
    :t₁ => t₁,
    :t₂ => t₂,
)

dt = 0.1 # s
model = StandardABM(Xie{3,4}, space, dt; properties)
````

A peculiarity of the `Xie` model is that the chemotactic properties of the
microbe differ between the forward and backward motile states, so we can
probe the response function in both the forward and backward motile state
by initializing two distinct microbes in the two states.
To keep the microbes in these motile states for the entire experiment duration,
we suppress their tumbles, and (just for total consistency with experiments)
we also set their speed to 0.

````@example 3_xie_response-function
add_agent!(model; motility=RunReverseFlick(Inf, [0], 0.0, [0]))
add_agent!(model; motility=RunReverseFlick(0.0, [0], Inf, [0]))
model[2].motility.current_state = 3 # manually set to backward run state
````

In addition to the `bias`, we will also monitor two other quantities
`state_m` and `state_z` which are internal variables of the `Xie` model
which represent the methylation and dephosphorylation processes which
together control the chemotactic response of the bacterium.

````@example 3_xie_response-function
nsteps = round(Int, T/dt)
adata = [bias, :state_m, :state_z]
adf, = run!(model, nsteps; adata)

S = Analysis.adf_to_matrix(adf, :bias)
m = (Analysis.adf_to_matrix(adf, :state_m))[:,1] # take only fw
z = (Analysis.adf_to_matrix(adf, :state_z))[:,1] # take only fw
````

We first look at the response function in the forward and backward
motile state: when the concentration increases we have a sharp negative
response (the tumble bias decreases), then the bacterium adapts to the new
concentration level, and when it drops back to the basal level we observe
a sharp positive response (the tumble bisa increases) before adapting
again to the new concentration level.

````@example 3_xie_response-function
_green = palette(:default)[3]
plot()
x = (0:dt:T) .- t₁
plot!(
    x, S,
    lw=1.5, lab=["Forward" "Backward"]
)
plot!(ylims=(-0.1,2.6), ylab="Response", xlab="time (s)")
plot!(twinx(),
    x, t -> concentration_field(t.+t₁,C₀,C₁,t₁,t₂),
    ls=:dash, lw=1.5, lc=_green, lab=false,
    tickfontcolor=_green,
    ylab="C (μM)", guidefontcolor=_green
)
````

By analysing the methylation and dephosphorylation processes, we can
understand how the chemotactic response arises.
First, when the concentration increases, both `m` and `z` increase and converge
to a new steady-state value, but since they respond on different timescales,
the response (defined by the difference between these two quantities), shows
a sharp decrease followed by a slower relaxation.
The same occurs for the negative stimulus.

````@example 3_xie_response-function
x = (0:dt:T) .- t₁
τ_m = model[1].adaptation_time_m
τ_z = model[1].adaptation_time_z
M = m ./ τ_m
Z = z ./ τ_z
R = M .- Z
plot(
    x, [M Z R],
    lw=2,
    lab=["m/τ_m" "z/τ_z" "m/τ_m - z/τ_z"],
    xlab="time (s)"
)
````

