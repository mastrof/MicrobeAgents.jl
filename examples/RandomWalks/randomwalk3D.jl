using MicrobeAgents
using Distributions
using Plots

## model parameters
L = 500
space = ContinuousSpace((L,L,L))
dt = 0.1
nsteps = 600

## abm setup
model = StandardABM(Microbe{3}, space, dt)
# add bacteria with different motile properties
add_agent!(model; motility=RunReverse(speed_forward=[55]), rotational_diffusivity=0.2)
add_agent!(model; motility=RunTumble(speed=Normal(30,6)), turn_rate=0.5)
add_agent!(model; motility=RunReverseFlick(speed_backward=[6]), rotational_diffusivity=0.1)

## simulation
adata = [:pos]
adf, _ = run!(model, nsteps; adata)

## postprocessing
traj = MicrobeAgents.unfold(vectorize_adf_measurement(adf,:pos), L)
x = first.(traj)
y = map(s -> s[2], traj)
z = last.(traj)
t = axes(x,1) .* dt

## visualization
plot(x, y, z, xlab="x", ylab="y", zlab="z", lw=2, ratio=1,
    lab=["RunReverse" "RunTumble" "RunReverseFlick"]
)
