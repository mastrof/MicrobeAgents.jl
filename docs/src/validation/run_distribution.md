```@meta
EditURL = "../../../examples/Validation/run_distribution.jl"
```

# Run-length distribution

In a random walk where the tumbling events occur as Poissonian events,
the resulting distribution of run lengths should be exponentially distributed,
with the average run length satisfying ``\tau = 1/\nu`` where ``\nu`` is
the average unbiased turn rate of the bacteria.

Since we don't work with continuous event-based simulations, the integration
timestep plays an important role in the resulting distribution.
If the timestep is too large compared to the typical run lengths,
the resulting distibution will be affected.
But if the timestep is too small, the simulations will be too expensive.
In most scenarios, a timestep ``\Delta t \sim \tau/10`` is typically the largest
value with which the correct run length distribution can be properly sampled.

````@example run_distribution
using MicrobeAgents

L = 1000
space = ContinuousSpace((L,L,L))
dt = 0.1 # s
model = StandardABM(Microbe{3}, space, dt)

runtime_expected = 2.0 # s
unbiased_turn_rate = 1 / runtime_expected # Hz
n = 500
for i in 1:n
    add_agent!(model; turn_rate=unbiased_turn_rate)
end

nsteps = 10000
adata = [velocity]
adf, = run!(model, nsteps; adata)
````

To verify that we sample the expected distribution, we will first
estimate the run lengths in our simulation with `Analysis.detect_turns!`
and `Analysis.run_durations`, and then we will fit the resulting distribution
with an `Exponential` distribution from Distributions.jl.
We expect the fitted distribution to match our histogram and return a
timescale `tau` close to the value of `runtime_expected`.

````@example run_distribution
using Distributions
Analysis.detect_turns!(adf)
run_lengths = vcat(Analysis.run_durations(adf)...) .* dt
estimated_pdf = fit(Exponential, run_lengths)
τ = scale(estimated_pdf)

using Plots
histogram(run_lengths; bins=dt/2:2dt:5τ, normalize=:pdf,
    lab="Simulation; τ = $(runtime_expected) s", lw=0)
plot!(dt:dt:5τ, t -> pdf(estimated_pdf, t-dt); lw=2,
    lab="exp(-t/τ)/τ; τ = $(round(τ;digits=3)) s")
plot!(xlab="Δt (s)", ylab="P(Δt)")
````

