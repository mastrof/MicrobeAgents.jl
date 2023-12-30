# # 1D Random walk

#=
Here we simulate a population of one-dimensional random walkers.

First we shall set up the model: we need to define the microbe type,
the space and the integration timestep.

For the microbe type, we will choose the base type `Microbe{1}`.

For the space, we need to set the domain size and whether it is periodic or not.
We will use a periodic box with an extent of 1000 μm;
if unspecified, `ContinuousSpace` will default to `periodic=true`.

For the integration timestep, we choose a value of 0.1 s.

Remember that lengths are always assumed to be in units of microns, times in seconds.
=#

using MicrobeAgents

L = 1000 # space size in μm
space = ContinuousSpace((L,))
dt = 0.1 # integration timestep in s
model = StandardABM(Microbe{1}, space, dt)

#=
Now that the model is initialized, we will add 10 microbes.
If we don't provide any argument on creation, default values from the constructor
will be used, i.e., an isotropic `RunTumble` motility, with speed 30 μm/s, an
unbiased tumbling rate of 1 Hz...

With the first (optional) argument to `add_agent!`, we can define the starting
position of the microbe. For convenience, we will initialize all of them
from position `(0,)`.
=#

n = 10 # number of microbes to add
foreach(_ -> add_agent!((0,), model), 1:n)

#=
We can now run the simulation.
We just need to define how many timesteps we want to simulate
and what kind of data we want to store during the simulation.
In this simulation, we only want to collect the microbe positions
at each timestep; we then set the `adata` vector (agent data)
to collect the position of the microbes.

The `run!` function will return a dataframe for the agent data (`adf`)
and one for the model data (`mdf`) collected during the simulation.
Here we are only collecting agent data and no model data.
=#

nsteps = 600
adata = [position]
adf, _ = run!(model, nsteps; adata);

#=
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
=#

Analysis.unfold!(adf, model)
traj = Analysis.adf_to_matrix(adf, :position_unfold)
x = first.(traj)
t = axes(x,1) .* dt

using Plots
plot(t, x, leg=false, xlab="time", ylab="displacement")
