# Random walks

Generating random walks with MicrobeAgents.jl is super-easy, and the walk properties
can be fine-tuned to match all sorts of needs.

By default, MicrobeAgents.jl produces random walks composed of
runs at constant speed and arbitrary reorientations where the waiting times
between such reorientations are i.i.d. random variables.
In the absence of chemotaxis (or other behavioral responses that affect microbe motility),
the generated random walks will display an exponential distribution of waiting times.

## Random Walk in D=1
```
using MicrobeAgents

L = 1000
space = ContinuousSpace((L,))
dt = 0.1
n = 10
nsteps = 600

model = StandardABM(Microbe{1}, space, dt; container=Vector)
foreach(_ -> add_agent!((0,), model), 1:n)
```
By default, `ContinuousSpace` will create a periodic domain; this will
allow us to mimic an "infinite" system.
All the microbes have been initialized from position 0, without
specifying any further property, so they will be initialized with a random
velocity (either `(+1,)` or `(-1,)` in this 1D scenario),
`speed=30.0`, and `turn_rate=1.0`.
The motility is set by default to `RunTumble(speed=[30.0])` but any other
motile pattern in 1D would produce the same result;
only speed is relevant here.

To run the simulation while collecting the bacterial positions
at each step we will then run
```
adata = [position]
adf, _ = run!(model, nsteps; adata)
```
Since the simulation box is periodic, the trajectories we collected in `adf`
fold around at the box edges. We can unfold them with the `unfold!` function
available through the `Analysis` submodule.
The unfolded coordinates will be stored in a new column `:position_unfold`.
For convenient plotting we can turn this new dataframe column into a matrix with
one trajectory per column with `adf_to_matrix`.
```
Analysis.unfold!(adf, model)
trajectories = Analysis.adf_to_matrix(adf, :position_unfold)
```
`trajectories` is now a `Matrix{SVector{1,Float64}}`.
To obtain the `x` positions for plotting we can call
```
x = first.(trajectories)
t = axes(x,1) .* dt
plot(t,x,lab=false,xlab="time",ylab="displacement")
```

## Random walks with different motile patterns in D=2
The procedure to generate a random walk in higher dimensions is
exactly the same
```
L = 500
space = ContinuousSpace((L,L))
dt = 0.1
nsteps = 600

model = StandardABM(Microbe{2}, space, dt; container=Vector)
```
But we can now add bacteria with different motility patterns
(we import Distributions.jl to use a Normal distribution for the speed of
one of our microbes)
```
using Distributions
add_agent!(model; motility=RunReverse(speed_forward=[55]), rotational_diffusivity=0.2)
add_agent!(model; motility=RunTumble(speed=Normal(30,6)), turn_rate=0.5)
add_agent!(model; motility=RunReverseFlick(speed_backward=[6]), rotational_diffusivity=0.1)
```
and then we run and visualize as before
```
adata = [:pos]
adf, _ = run!(model, nsteps; adata)

# postprocessing
traj = MicrobeAgents.unfold(vectorize_adf_measurement(adf,:pos), L)
x = first.(traj)
y = last.(traj)
t = axes(x,1) .* dt
plot(x, y, xlab="x", ylab="y", ratio=1,
    lab=["RunReverse" "RunTumble" "RunReverseFlick"]
)
```
![Two-dimensional random walks with different motility patterns](rw2d.svg)
