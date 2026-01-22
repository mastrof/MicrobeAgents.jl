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
3. Initialize an `AgentBasedModel` object with the desired space, microbe type, integration timestep, and any extra property needed for the simulation.
4. Populate the ABM with microbe instances.
5. Choose the observables to collect during production and run the model.

MicrobeAgents.jl re-exports and extends various function from Agents.jl in order to work
as a standalone, but it is generally recommended to use it in combination with
Agents.jl for extra goodies.

## Space
MicrobeAgents.jl only supports continuous spaces with dimensions 1, 2 or 3.
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
Microbes are represented by subtypes of the `AbstractMicrobe` type, which is in turn a subtype of `AbstractAgent` from Agents.jl
```@docs
AbstractMicrobe
```

MicrobeAgents provides different `AbstractMicrobe` subtypes representing different models of bacterial behavior from the scientific literature.
The list of implemented models can be obtained with `subtypes(AbstractMicrobe)`.

A basic type, which is typically sufficient for simple motility simulations and does not include chemotaxis, is the `Microbe` type.
```@docs
Microbe
```
Microbe instances should only be created within an `AgentBasedModel`, the fundamental structure which embeds everything that has to do with the agent-based simulations you want to run.
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

Now, with the `add_agent!` function we will populate the model with microbes of
the specified type (`Microbe{1}`).
The only argument we must *always* specify for `add_agent!` is the motility of the microbe, via the `motility` keyword. An overview of the motility interface is given later;
For now we will just use a Run-Tumble motility with an average run duration of 1 second and a constant swimming speed of 20 micron / second.
The third required argument to the call is the distribution of reorientation
angles, but it is irrelevant in 1-dimensional systems so we can just pass
an empty array.
If unspecified, position, direction and speed of the microbe will be assigned randomly;
all the other fields will be assigned default values from the constructor (unless specified).
To select a position, it can be passed as the first argument to the `add_agent!` call,
and any other bacterial parameter can be defined via keyword arguments.
All of the following are valid calls
```
motility = RunTumble(1.0, [20.0], [])
# a Microbe with large radius
add_agent!(model; motility, radius=10.0)
# a Microbe with custom position and high coefficient of rotational diffusion
add_agent!((53.2,), model; motility, rotational_diffusivity=0.5)
# a Microbe initialized with velocity to the right
add_agent!(model; motility, vel=(1.0,))
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

The `microbe_step!` function is split into 3
subroutines, corresponding to movement, internal state changes, and active reorientations, applied in this order:

```@docs
move_step!
affect_step!
reorient_step!
```

All of them are exported and may be overwritten with user-specified
behavior, or their order within `microbe_step!` changed
(e.g. to perform internal state updates before movement),
which is generally more convenient than fully rewriting
`microbe_step!` in case modifications need to be applied.


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
Since `adf.position` is a vector of static vectors,
we first have to unpack the x and y values
and then we are ready to plot our trajectory.
```
using Plots
x = first.(adf.position)
y = last.(adf.position)
plot(x, y)
```


## Motility patterns
In MicrobeAgents.jl, motility patterns are represented through the
`Motility{N}` type.
A `Motility` is composed of `N` instances of `MotileState` and of a set
of transition probabilities between these `N` states.
A `MotileState` contains information about the speed distribution,
turn angle distributions, and average lifetime of a particular
state of motion.

There are two kinds of `MotileState`s, `RunState` and `TurnState`.
A `RunState` is used to represent states associated with translational
motion and no angular reorientation (e.g. `Run`).
A `TurnState` is instead used to represent states where no translational
motion happens, but where reorientations may occur
(e.g. `Tumble`, `Reverse`, `Flick`, `Stop`).

Instances of `MotileState`s can be arbitrarily combined into a `Motility`.
```@docs
Motility
```
The most common motility patterns are, however, pre-implemented.
```@docs
RunTumble
RunReverse
RunReverseFlick
RunStop
```
