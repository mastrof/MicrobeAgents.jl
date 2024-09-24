```@meta
EditURL = "../../../../examples/RandomWalks/1_randomwalk1D.jl"
```

# 1D Random walk

Here we simulate a population of one-dimensional random walkers.

First we shall set up the model: we need to define the microbe type,
the space and the integration timestep.

For the microbe type, we will choose the base type `Microbe{1}`, where
parameter `1` refers to the number of dimensions in which the microbe can move.

For the space, we need to set the domain size and whether it is periodic or not.
We will use a periodic box with an extent of 1000 μm;
if unspecified, `ContinuousSpace` will default to `periodic=true`.

For the integration timestep, we choose a value of 0.1 s.

Remember that lengths are always assumed to be in units of microns, times in seconds.

````@example 1_randomwalk1D
using MicrobeAgents

L = 1000 # space size in μm
space = ContinuousSpace((L,))
dt = 0.1 # integration timestep in s
model = StandardABM(Microbe{1}, space, dt)
````

Now that the model is initialized, we will add 10 microbes.

With the first (optional) argument to `add_agent!`, we can define the starting
position of the microbe. For convenience, we will initialize all of them
from position `(0,)`.

Moreover, it's always required to specify the motility of the microbes via
the keyword argument `motility`.
For example, `RunTumble` needs as inputs the average duration of runs,
the distribution of run speeds and the distribution of reorientation angles.
The latter is irrelevant in 1D, since the bacterium can only revert its
direction along the line, but we still have to pass it.
We can just use an empty vector.
Passing a one-element vector `[U]` as the speed distribution
means that all runs will have the same velocity `U`.
We could have also specified the average duration of tumbles
(equal to reversals in 1 dimension) via an optional 4th argument (try it! e.g. 0.5);
when unspecified, the tumbles are taken to be instantaneous.

````@example 1_randomwalk1D
n = 10 # number of microbes to add
τ_run = 1.0 # average run duration in s
U = 30.0 # swimming speed in μm/s
motility = RunTumble(τ_run, [U], [], 0.5)
foreach(_ -> add_agent!((0,), model; motility), 1:n)
````

We can now run the simulation.
We just need to define how many timesteps we want to simulate
and what kind of data we want to store during the simulation.
In this simulation, we only want to collect the microbe positions
at each timestep; we then set the `adata` vector (agent data)
to collect the position of the microbes.

The `run!` function will return a dataframe for the agent data (`adf`)
and one for the model data (`mdf`) collected during the simulation.
Here we are only collecting agent data and no model data.

````@example 1_randomwalk1D
nsteps = 600
adata = [position]
adf, _ = run!(model, nsteps; adata);
nothing #hide
````

The simulation is done. We now want to visualize our results.
One last thing we need to is to "unfold" the trajectories of our microbes.
In fact, since we used a periodic domain, if we just plotted the trajectories
we would see them crossing between the two sides of the box, which is not what we want.

With the unfolding, the trajectories are expanded as if they were simulated in an
infinite system.

The `unfold!` function and other useful post-processing functions can be
accessed through the `Analysis` submodule of MicrobeAgents.jl.
`unfold!` operates directly on the agent dataframe `adf`
and generates a new column `:position_unfold` from the raw `:position` column.
After the unfolding, we extract the individual trajectories in the form of a
matrix with `adf_to_matrix`, where each trajectory will be stored as a column.
This will be convenient for plotting.

````@example 1_randomwalk1D
Analysis.unfold!(adf, model)
traj = Analysis.adf_to_matrix(adf, :position_unfold)
x = first.(traj)
t = axes(x,1) .* dt

using Plots
plot(t, x, leg=false, xlab="time", ylab="displacement")
````

