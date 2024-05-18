```@meta
EditURL = "../../../examples/Validation/velocity_autocorrelations.jl"
```

# Velocity autocorrelation functions

Orientational correlations in the microbe's motion can be
probed through velocity autocorrelation functions.
Each motility pattern is characterized by a specific form
of the velocity autocorrelation function, which have been analytically
evaluated by Taktikos *et al* (PLoS ONE 2012).
We will compare the velocity autocorrelation functions of
run-tumble, run-reverse and run-reverse-flick motility, in the absence
of rotational diffusion, to the analytical calculations:
```math
\phi(t) = \begin{cases}
    \exp(-t/\tau), & \text{run-tumble} \\
    \exp(-2t/\tau), & \text{run-reverse} \\
    (1 - t/2\tau)\exp(-t/\tau), & \text{run-reverse-flick}
\end{cases}
```

````@example velocity_autocorrelations
using MicrobeAgents
using Plots

U = 30.0 # μm/s
τ_run = 1.0 # s
turn_rate = 1 / τ_run # 1/s
Δt = 0.01 # s
L = 1e4 # μm
space = ContinuousSpace((L,L,L))

model = StandardABM(Microbe{3}, space, Δt; container=Vector)
n = 200
for Mot in (RunTumble, RunReverse, RunReverseFlick), i in 1:n
    if Mot == RunTumble
        motility = RunTumble(τ_run, [U], 0.0)
    else
        motility = Mot(τ_run, [U], τ_run, [U])
    end
    add_agent!(model; motility)
end
# keep track of ids of each motile strategy
ids_runtumble = 1:n
ids_runreverse = (1:n) .+ n
ids_runrevflick = (1:n) .+ 2n

nsteps = round(Int, 100τ_run / Δt)
adata = [direction]
adf, = run!(model, nsteps; adata)

# separate the dataframes by motile strategy
adf_runtumble = filter(:id => id -> id in ids_runtumble, adf; view=true)
adf_runrev = filter(:id => id -> id in ids_runreverse, adf; view=true)
adf_runrevflick = filter(:id => id -> id in ids_runrevflick, adf; view=true)
adfs = [adf_runtumble, adf_runrev, adf_runrevflick]

# evaluate the autocorrelation functions
using StatsBase: mean
t = range(0, (nsteps-1)*Δt; step=Δt)
Φ = hcat([mean(Analysis.acf(a, :direction; normalize=true)) for a in adfs]...)

# from Taktikos et al.
Φ_theoretical = hcat([
    exp.(-t ./ τ_run),
    exp.(-t ./ (τ_run / 2)),
    (1 .- t ./ (2τ_run)) .* exp.(-t ./ τ_run),
]...)

plot(
    xlims=(0,6τ_run), ylims=(-0.1, 1.05),
    xlab="Δt / τ",
    ylab="velocity autocorrelation",
)
plot!(t, Φ_theoretical, lw=2, lc=[1 2 3],
    label=["Run-Tumble" "Run-Reverse" "Run-Reverse-Flick"]
)
# show simulation values only at selected lags for better visualization
lags = range(0, length(t)-1; step=20)
scatter!(t[lags.+1], Φ[lags.+1,:], m=:x, mc=[1 2 3], msw=2, label=false)
hline!([0.0], lw=0.8, ls=:dash, lc=:black, lab=false)
````

