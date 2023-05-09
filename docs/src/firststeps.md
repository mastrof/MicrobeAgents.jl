# First steps
## Creating a microbe
Microbes are represented by subtypes of the `AbstractMicrobe` type, which is in turn a subtype of `AbstractAgent` introduced by Agents.jl
```@docs
AbstractMicrobe
```

MicrobeAgents provides different `AbstractMicrobe` subtypes representing different models of bacterial behavior from the literature.

A basic type, which is typically sufficient for simple motility simulations and does not include chemotaxis, is the `Microbe` type.
```@docs
Microbe
```

The dimensionality of `Microbe` *must* always be specified on creation. All the fields are instead optional, and if not specified will be assigned default values.

To create a `Microbe` living in a 1-dimensional space, with run-tumble motility and average turn rate ``\nu=1\;s^{-1}``, it is therefore sufficient to run
```
Microbe{1}()
```
Similarly, for two and three dimensions:
```
Microbe{2}()
Microbe{3}()
```

Any custom parameter can be set via kwargs:
```
Microbe{3}(
    turn_rate = 0.6,
    rotational_diffusivity = 0.1
)
```


All the other subtypes of `AbstractMicrobe` work in a similar way, although
they will have distinct default values and extra fields.

```@docs
BrownBerg
Brumley
Celani
Xie
```


## Creating a model
MicrobeAgents.jl exploits the `AgentBasedModel` interface from Agents.jl.
While the standard Agents.jl syntax will always work, it is typically more
convenient to use the method extensions provided by MicrobeAgents.jl, which
also includes some default parameters required by the simulations.
Both `StandardABM` and `UnremovableABM` are supported.
Whenever removal of microbes during the simulation is not needed,
`UnremovableABM` is the recommended choice.
```@docs
UnremovableABM
StandardABM
```

To create a simple model, we just need to choose a microbe type, the size of
the simulation domain and the integration timestep.
The properties of the simulation domain are wrapped in the `ContinuousSpace`
object (re-exported from Agents.jl).
```
extent = (1000.0, 500.0) # size of 2D simulation domain
space = ContinuousSpace(extent)
dt = 0.1 # integration timestep
model = UnremovableABM(Microbe{2}, space, dt)
```
By default this creates a model with periodic boundary conditions.
For hard wall boundary conditions we can instead specify
`space=ContinuousSpace(extent; periodic=false)`.

Now bacteria can be added with `add_agent!` function.
```@docs
add_agent!
```

We need not specify anything if we want the microbe to be added at a random
position with the default values from the constructor.
```
add_agent!(model)
```
The microbe will be now accessible as `model[1]`.

The ABM is now ready to run. We just need to specify how many steps we want
to simulate and what data to collect during the run:
```
nsteps = 100
adf, mdf = run!(model, nsteps; adata=[:pos])
```
`run!` will return two dataframes, one for the agent-level data (`adf`) and one
for the model-level data (`mdf`, which in this case will be empty).
This way, we have produced our first random walk.
Since `adf.pos` is a vector of tuples, we first have to unpack the x and y values
and then we are ready to plot our trajectory.
```
using Plots
x = first.(adf.pos)
y = last.(adf.pos)
plot(x, y)
```

## Motility patterns
In MicrobeAgents.jl, motility patterns are represented as instances of
`AbstractMotility`.
In particular, currently available patterns are distinguished into two further
categories: `AbstractMotilityOneStep` or `AbstractMotilityTwoStep`.
```@docs
AbstractMotilityOneStep
AbstractMotilityTwoStep
```
One-step motility pattern are characterized by a single swimming stage. Two-step
motility patterns instead have two stages which can have distinct properties;
these two stages are referred to as "forward" and "backward" but can really
represent anything.

MicrobeAgents.jl defines three standard motility patterns:
```@docs
RunTumble()
RunReverse()
RunReverseFlick()
```
The default values provided by the constructors always consider a constant
swimming speed of 30 (micron/s) and "ideal" turn-angle distributions.
For more accurate simulation where the reorientation statistics of the microbes
is important, appropriate distributions should be specified;
the constructors will accept any object that can be sampled via `rand()`.
