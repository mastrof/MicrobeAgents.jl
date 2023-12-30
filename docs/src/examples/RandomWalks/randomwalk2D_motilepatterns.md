```@meta
EditURL = "../../../../examples/RandomWalks/randomwalk2D_motilepatterns.jl"
```

# 2D Random walk and motile patterns

Here we will simulate two dimensional random walk with different motile patterns.

As usual we start by setting up the model.

````@example randomwalk2D_motilepatterns
using MicrobeAgents
using Distributions
using Plots

L = 500 # space size in μm
space = ContinuousSpace((L,L)) # defaults to periodic
dt = 0.1 # integration timestep in s
model = StandardABM(Microbe{2}, space, dt)
````

We will now add microbes individually, choosing different properties for each.
The motile pattern can be customized through the `motility` keyword.

The `RunTumble` motility consists of straight runs interspersed with isotropic
reorientations (tumbles). We can define the `speed` of follow a a `Normal` distribution
(from Distributions.jl) with mean 30 μm/s and standard deviation 6 μm/s.
This means that, after every tumble, the microbe will change its speed following
this distribution.
Further, we reduce the unbiased tumbling rate (`turn_rate`) of the microbe
from the default value of 1 Hz to 0.5 Hz.

````@example randomwalk2D_motilepatterns
add_agent!(model; motility=RunTumble(speed=Normal(30,6)), turn_rate=0.5)
````

The `RunReverse` motility consists of alternating straight runs and 180-degree reversals.
Differently from the `RunTumble`, the `RunReverse` can be considered as a two-step motility
pattern and we can assign different properties to the "forward" and the "backward" state of motion.
If no properties are explicitly specified for the "backward" state, it will inherit those of
the "forward" state.
Here we just set the speed to the constant value of 55 μm/s.
Further, we set the `rotational_diffusivity` of the microbe to 0.2 rad²/s; in absence of
rotational diffusion, the run reverse motility is pathologically incapable of exploring space efficiently.

````@example randomwalk2D_motilepatterns
add_agent!(model; motility=RunReverse(speed=[55]), rotational_diffusivity=0.2)
````

The `RunReverseFlick` motility consists of a straight run, a 180-degree reversal, then another
straight run followed by a 90-degree reorientation (the flick).
Like the `RunReverse`, this is a two-step motility pattern; indeed, we can imagine it as
a `RunReverse` where the reorientation in the "backward" motile state is 90 instead of 180 degrees.
We set the `speed_backward` to 6 μm/s, while the speed in the forward mode will keep its default
value (30 μm/s). We also set the rotational diffusivity to 0.1 rad²/s.

````@example randomwalk2D_motilepatterns
add_agent!(model; motility=RunReverseFlick(speed_backward=[6]), rotational_diffusivity=0.1)
````

Now we can run (collecting the microbe positions at each timestep), unfold the trajectories,
and visualize them.

````@example randomwalk2D_motilepatterns
nsteps = 600
adata = [position]
adf, _ = run!(model, nsteps; adata)

Analysis.unfold!(adf, model)
traj = Analysis.adf_to_matrix(adf, :position_unfold)
x = first.(traj)
y = last.(traj)
t = axes(x,1) .* dt

plot(x, y, xlab="x", ylab="y", ratio=1,
    lab=["RunTumble" "RunReverse" "RunReverseFlick"]
)
````

