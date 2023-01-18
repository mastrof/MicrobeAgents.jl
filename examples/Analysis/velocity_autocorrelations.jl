using MicrobeAgents
using Plots

U = 30.0 # μm/s
τ_run = 1.0 # s
turn_rate = 1 / τ_run # 1/s
Δt = 0.01 # s
L = 1e4 # μm
extent = (L,L,L)

model = ABM(Microbe{3}, extent, Δt)
n = 200
for i in 1:n
    add_agent!(model; turn_rate, motility=RunTumble(speed=[U]))
    add_agent!(model; turn_rate, motility=RunReverse(speed_forward=[U]))
    add_agent!(model; turn_rate, motility=RunReverseFlick(speed_forward=[U]))
end

nsteps = round(Int, 100τ_run / Δt)
adata = [:vel]
adf, = run!(model, nsteps; adata)

adf_runtumble = filter(:id => id -> model.agents[id].motility isa RunTumble, adf; view=true)
adf_runrev = filter(:id => id -> model.agents[id].motility isa RunReverse, adf; view=true)
adf_runrevflick = filter(:id => id -> model.agents[id].motility isa RunReverseFlick, adf; view=true)
adfs = [adf_runtumble, adf_runrev, adf_runrevflick]

lags = unique(round.(Int, range(0, findfirst(t.>6), length=30)))
@time Φ = hcat([acf(a,:vel,lags) for a in adfs]...)

t = range(0, (nsteps-1)*Δt; step=Δt)
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
scatter!(t[lags.+1], Φ ./ U^2, m=:x, mc=[1 2 3], label=false)
hline!([0.0], lw=0.8, ls=:dash, lc=:black, lab=false)