# First steps

## Structure
An `AgentBasedModel` object embeds all the properties of the system to be simulated
and maps unique IDs to microbe instances.
During the simulation, the model is evolved in discrete time steps, with
each microbe's "state" being updated according to specified rules.
Standard rules for motion, reorientations and chemotaxis are available by default,
but custom behaviors can be implemented via user-defined functions.

The typical workflow to run a simulation in MicrobeAgents.jl goes as follows:
1. Define the size and properties of the space in which the microbes will move.
2. Choose an appropriate microbe type to represent the desired behavior, or define a new one.
3. Initialize an `AgentBasedModel` object with the desired space, microbe type, integration time step, and any extra property needed for the simulation.
4. Populated the ABM with microbe instances.
5. Run the model (defining custom stepping functions if required) and collect data.

MicrobeAgents.jl re-exports and extends various function from Agents.jl in order to work
as a standalone, but it is generally recommended to use it in combination with
Agents.jl for extra goodies.

## Space
MicrobAgents.jl only supports continuous spaces with dimensions 1, 2 or 3.
Spaces can be created with the `ContinuousSpace` function (reexported from Agents.jl).
The extent of the space must be given as a tuple, and periodicity is set with
the `periodic` kwarg (defaults to true).
```
# one-dimensional periodic space
extent = (1.0,)
ContinuousSpace(extent)

# two-dimensional non-periodic space
extent = (1.0, 2.0)
ContinuousSpace(extent; periodic=false)
```


## Microbes
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

To create a `Microbe` living in a 1-dimensional space, with default parameters
(`RunTumble` motility, average turn rate ``\nu=1\;s^{-1}``, and no rotational diffusivity),
it is therefore sufficient to run
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
Default values are typically assigned following the original implementation in the literature.

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
Whenever removal of microbes during the simulation is not needed,
it is recommended to call `StandardABM` with the `container=Vector`
keyword argument to improve performance.
```@docs
StandardABM
```

To create a simple model, we just need to choose a microbe type, the size of
the simulation domain and the integration timestep.
The properties of the simulation domain are wrapped in the `ContinuousSpace`
object.
```
extent = (1000.0, 500.0) # size of 2D simulation domain
space = ContinuousSpace(extent)
dt = 0.1 # integration timestep
model = StandardABM(Microbe{2}, space, dt)
```

Now bacteria can be added with the `add_agent!` function.
```@docs
add_agent!
```

We need not specify anything if we want the microbe to be added at a random
position with the default values from the constructor.
```
add_agent!(model)
```
The microbe will be now accessible as `model[1]`.

## Running a model
After the model has been created and populated with the desired number of microbes,
we are ready to run the simulation.
We just need to specify how many steps we want
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

Notice that we did not specify at all how the timestepping is performed.
MicrobAgents.jl implements a default timestepper which is applied to all
`AbstractMicrobe` instances, which takes care of motion, rotational diffusion
and reorientations.
Each subtype is then equipped with its own `affect!` and `turnrate` functions
(explained later) which determine extra behavioral features (such as chemotaxis).

If different behavior is desired, the integration can be customized by
passing custom timestepping functions to `run!`:
```
run!(model, my_agent_step!, my_model_step!, nsteps)
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
swimming speed of 30 (micron/s) and "ideal" turn-angle distributions
(isotropic for tumbles, perfect 180 and 90 degree reorientations for reverse and flicks).
For more accurate simulation where the reorientation statistics of the microbes
is important, appropriate distributions should be specified;
the constructors will accept any object that can be sampled via `rand()`.
