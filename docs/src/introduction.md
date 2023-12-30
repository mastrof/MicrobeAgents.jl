# Introduction

## Structure
An `AgentBasedModel` object embeds all the properties of the system to be simulated
and maps unique IDs to microbe instances.
During the simulation, the model is evolved in discrete time steps, with
each microbe's position, velocity and "state" being updated according to specified rules.
Standard rules for motion, reorientations and chemotaxis are available by default,
but custom behaviors can be implemented via user-defined functions.

The typical workflow to run a simulation in MicrobeAgents.jl goes as follows:
1. Define the size and properties of the space in which the microbes will move.
2. Choose an appropriate microbe type to represent the desired behavior, or define a new one.
3. Initialize an `AgentBasedModel` object with the desired space, microbe type, integration time step, and any extra property needed for the simulation.
4. Populated the ABM with microbe instances.
5. Choose the observables to collect during production and run the model.

MicrobeAgents.jl re-exports and extends various function from Agents.jl in order to work
as a standalone, but it is generally recommended to use it in combination with
Agents.jl for extra goodies.

## Space
MicrobAgents.jl only supports continuous spaces with dimensions 1, 2 or 3.
Spaces can be created with the `ContinuousSpace` function (reexported from Agents.jl).
The extent of the space must be given as a tuple or `SVector`, and periodicity is set with
the `periodic` kwarg (defaults to true).
```
# one-dimensional periodic space
extent = (1.0,)
ContinuousSpace(extent)

# two-dimensional non-periodic space
extent = (1.0, 2.0)
ContinuousSpace(extent; periodic=false)

# three-dimensional space periodic only along the x direction
extent = (100.0, 20.0, 20.0)
ContinuousSpace(extent; periodic=(true,false,false))
```


## Microbes
Microbes are represented by subtypes of the `AbstractMicrobe` type, which is in turn a subtype of `AbstractAgent` introduced by Agents.jl
```@docs
AbstractMicrobe
```

MicrobeAgents provides different `AbstractMicrobe` subtypes representing different models of bacterial behavior from the literature.
The list of implemented models can be obtained with `subtypes(AbstractMicrobe)`.

A basic type, which is typically sufficient for simple motility simulations and does not include chemotaxis, is the `Microbe` type.
```@docs
Microbe
```

The dimensionality of `Microbe` *must* always be specified on creation. All the fields are instead optional, and if not specified will be assigned default values.
Microbe instances should only be created within an `AgentBasedModel`.

In MicrobeAgents.jl, models are created through the `StandardABM` function.
```@docs
StandardABM
```
To initialize a model we must provide the microbe type, the simulation space, and the
integration timestep of the simulation. All other parameters are optional.
To setup a model for `Microbe`s living in a 1-dimensional space we can therefore run
```
space = ContinuousSpace((100.0,); periodic=false)
dt = 0.1
model = StandardABM(Microbe{1}, space, dt)
```

Now, calling `add_agent!(model)` will populate the model with microbes of
the specified type (`Microbe{1}`) using the default values of the constructor,
and automatically generating a random position and a random velocity vector.
To select a position, it can be passed as the first argument to the `add_agent!` call,
and any other bacterial parameter can be defined via keyword arguments.
All of the following are valid calls
```
# a Microbe with large radius and low tumble rate
add_agent!(model; radius=10.0, turn_rate=0.17)
# a Microbe with custom position and high coefficient of rotational diffusion
add_agent!((53.2,), model; rotational_diffusivity=0.5)
# a Microbe initialized with velocity to the right
add_agent!(model; vel=(1.0,))
```

All the other subtypes of `AbstractMicrobe` work in a similar way, although
they will have distinct default values and extra fields.
When possible, default values are typically assigned following the original implementation in the literature.

```@docs
BrownBerg
Brumley
Celani
Xie
```


## More about models
MicrobeAgents.jl exploits the `AgentBasedModel` interface from Agents.jl.
While the standard Agents.jl syntax will always work, it is typically more
convenient to use the method extensions provided by MicrobeAgents.jl, which
also includes some default parameters required by the simulations.
If the simulation requires removal/addition of microbes, it is recommended
to call `StandardABM` with the `container=Dict` keyword argument,
otherwise MicrobeAgents.jl defaults to `container=Vector` which provides
better performance.

In addition to the microbe instances, the model should also wrap all
the other information required to perform the simulation.

MicrobeAgents.jl defines default timestepping functions which are used
to evolve the microbes and the model, and are accessible through the
`microbe_step!` and `model_step!` keywords in `StandardABM`.
By default, the `microbe_step!` function performs, in order:
- update microbe position according to current velocity
- randomize the microbe orientation through rotational diffusion (if present)
- update internal state of the microbe (e.g. chemotaxis or other user-defined behavior)
- perform reorientation events following Poissonian statistics
The `model_step!` function instead defaults to a dummy function which does nothing.
Any custom behavior can be implemented by simply modifying these two functions.

Any type of external parameter that should be used during the simulation should be
passed to `StandardABM` through the `properties` dictionary.

## Running a model
After the model has been created and populated with the desired number of microbes,
we are ready to run the simulation.
We just need to specify how many steps we want
to simulate and what data to collect during the run:
```
nsteps = 100
adf, mdf = run!(model, nsteps; adata=[position])
```
`run!` will return two dataframes, one for the agent-level data (`adf`) and one
for the model-level data (`mdf`, which in this case will be empty).
This way, we have produced our first random walk.
Since `adf.position` is a vector of tuples, we first have to unpack the x and y values
and then we are ready to plot our trajectory.
```
using Plots
x = first.(adf.position)
y = last.(adf.position)
plot(x, y)
```


## Motility patterns
In MicrobeAgents.jl, motility patterns are represented as instances of
`AbstractMotility`.
In particular, currently available patterns are distinguished into two further
categories: `AbstractMotilityOneStep` or `AbstractMotilityTwoStep`.
```@docs; canonical=false
MotilityOneStep
MotilityTwoStep
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
