using MicrobeAgents
using Plots

## model parameters
L = 1000
space = ContinuousSpace((L,))
dt = 0.1
n = 10
nsteps = 600

## abm setup
model = UnremovableABM(Microbe{1}, space, dt)
foreach(_ -> add_agent!((0,), model), 1:n)

## simulation
adata = [:pos]
adf, _ = run!(model, nsteps; adata)

## postprocessing
traj = MicrobeAgents.unfold(vectorize_adf_measurement(adf,:pos), L)
x = first.(traj)
t = axes(x,1) .* dt

## visualization
plot(t, x, leg=false, xlab="time", ylab="displacement")
