# # Run durations

#=
Switching between motile states is a Poisson process: a microbe in a run
state of average duration ``\tau`` leaves it with probability ``\Delta t/\tau``
at each timestep. For ``\Delta t \ll \tau`` the run durations are therefore
exponentially distributed:
```math
P(t) = \frac{1}{\tau} \exp(-t/\tau)
```
Behaviors can bias the switching rate by a multiplicative factor ``b``
(see `bias`), which rescales the mean run duration to ``\tau/b``.
We will compare the run-length distributions of unbiased microbes and of
microbes with a constant bias against the exponential prediction.
=#

using MicrobeAgents
using Plots

U = 30.0 # μm/s
τ_run = 1.0 # s
Δt = 0.01 # s
L = 1e4 # μm
n = 200 # microbes per population
nsteps = round(Int, 100τ_run / Δt)

space = ContinuousSpace((L, L, L))
model = StandardABM(Microbe{3}, space, Δt; container=Vector)
## instantaneous isotropic tumbles, no rotational diffusion:
## every change of direction marks the end of a run
motility = RunTumble([U], τ_run, Isotropic(3))
biases = (1.0, 2.0)
for b in biases, _ in 1:n
    add_agent!(model; motility, behaviors=(Behavior(bias = (m, model) -> b),))
end
adf, = run!(model, nsteps; adata=[velocity])
Analysis.detect_turns!(adf)
durations = Analysis.run_durations(adf)

plt = plot(
    xlab="run duration (s)", ylab="probability density",
    xlims=(0, 4τ_run),
)
for (k, b) in enumerate(biases)
    ids = (1:n) .+ (k-1)*n
    ## the first run of each microbe started before the observation: discard it
    T = Δt .* reduce(vcat, [d[2:end] for d in durations[ids]])
    τ = τ_run / b
    stephist!(plt, T, normalize=:pdf, bins=0:10Δt:10τ_run, lw=2, lc=k, label="bias = $b")
    plot!(plt, t -> exp(-t/τ)/τ, 0, 4τ_run, lw=2, ls=:dash, lc=:black, label=(k == 1 ? "theory" : false))
end
plt
