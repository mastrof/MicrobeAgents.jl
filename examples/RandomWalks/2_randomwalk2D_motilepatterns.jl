# # 2D Random walk and motile patterns

#=
Here we will simulate two dimensional random walk with different motile patterns.

As usual we start by setting up the model.
=#

using MicrobeAgents
using Distributions
using Plots

L = 500 # space size in μm
space = ContinuousSpace((L,L)) # defaults to periodic
dt = 0.1 # integration timestep in s
model = StandardABM(Microbe{2}, space, dt)

#=
We will now add microbes individually, choosing different properties for each.
The motile pattern can be customized through the `motility` keyword.

The `RunTumble` motility consists of straight runs interspersed with
reorientations (tumbles).
With the first argument, we set the average duration of tumbles, 1 s in this case.
We can then define the speed to follow a a `Normal` distribution
(from Distributions.jl) with mean 30 μm/s and standard deviation 6 μm/s.
This means that, after every tumble, the microbe will change its speed following
this distribution.
We should then specify the distribution of angular reorientations.
For convenience, MicrobeAgents implements an `Isotropic` function
which produces the appropriate distribution to have isotropic reorientations
in a space with given dimensionality.
In this case, using `Isotropic(2)` produces a `Uniform` distribution
of angles between -π and +π.
Finally, we can set tumbles to also have a finite duration, let's say 0.1 s.
=#
add_agent!(model; motility=RunTumble(1.0, Normal(30,6), Isotropic(2), 0.1))

#=
The `RunReverse` motility consists of alternating straight runs and 180-degree reversals.
This motility pattern is implemented as a 4-step motility (i.e., `Motility{4}`):
a run "forward", a reversal of direction, a run "backward",
and another reversal of direction, after which the cycle starts again.
In principle, each of these 4 steps can be modified independently from the others.
Here, we will set the duration and speed of the two run steps, and keep default
values for reversals.
Further, we set the `rotational_diffusivity` of the microbe to 0.2 rad²/s; in absence of
rotational diffusion, the run reverse motility is pathologically incapable of
exploring space efficiently.
=#
add_agent!(model; motility=RunReverse(1.0, [30.0], 0.7, [20.0]), rotational_diffusivity=0.2)

#=
The `RunReverseFlick` motility consists of a straight run, a 180-degree reversal, then another
straight run followed by a 90-degree reorientation (the flick).
Like the `RunReverse`, this is a four-step motility pattern; indeed, we can imagine it as
a `RunReverse` where the reorientation after the "backward" run is 90 instead of 180 degrees.
Again, we will set only run durations and speeds.
We also set the rotational diffusivity to 0.1 rad²/s.
=#
add_agent!(model; motility=RunReverseFlick(2.0, [25.0], 0.5, [25.0]), rotational_diffusivity=0.1)

#=
Now we can run (collecting the microbe positions at each timestep), unfold the trajectories,
and visualize them.
=#

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
