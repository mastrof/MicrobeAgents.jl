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
for Motility in (RunTumble, RunReverse, RunReverseFlick), i in 1:n
    add_agent!(model; turn_rate, motility=Motility(speed=[U]))
end
ids_runtumble = 1:n
ids_runreverse = (1:n) .+ n
ids_runrevflick = (1:n) .+ 2n

nsteps = round(Int, 100τ_run / Δt)
adata = [:vel]
adf, = run!(model, nsteps; adata)

adf_runtumble = filter(:id => id -> id in ids_runtumble, adf; view=true)
adf_runrev = filter(:id => id -> id in ids_runreverse, adf; view=true)
adf_runrevflick = filter(:id => id -> id in ids_runrevflick, adf; view=true)
adfs = [adf_runtumble, adf_runrev, adf_runrevflick]

t = range(0, (nsteps-1)*Δt; step=Δt)
lags = unique(round.(Int, range(0, findfirst(t.>6), length=30)))
@time Φ = hcat([acf(a,:vel,lags) for a in adfs]...)

Φ_theoretical = hcat([
    exp.(-t ./ τ_run),
    exp.(-t ./ (τ_run / 2)),
    (1 .- t ./ (2τ_run)) .* exp.(-t ./ τ_run),
]...) # Taktikos et al. 2013 PLoS ONE

plot(
    xlims=(0,6τ_run), ylims=(-0.1, 1.05),
    xlab="Δt / τ",
    ylab="velocity autocorrelation",
)
plot!(t, Φ_theoretical, lw=2, lc=[1 2 3], label=["Run-Tumble" "Run-Reverse" "Run-Reverse-Flick"])
scatter!(t[lags.+1], Φ, m=:x, mc=[1 2 3], msw=2, label=false)
hline!([0.0], lw=0.8, ls=:dash, lc=:black, lab=false)
