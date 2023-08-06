using MicrobeAgents
using Distributions
using Plots

L = 1000
space = ContinuousSpace((L,L))
dt = 0.01

model = UnremovableABM(Microbe{2}, space, dt)
n = 500
Drot = 0.1
τ = 1.5
for i in 1:n
    add_agent!(model; turn_rate=1/τ, rotational_diffusivity=Drot)
end

nsteps = 10000
adf, _ = run!(model, nsteps; adata=[:vel])

ᾱ = 4 * sqrt(2*Drot*dt)
runs = rundurations(adf, dt; threshold_angle=ᾱ)
histogram(vcat(runs...), bins=dt/2:10dt:5τ, normalize=:pdf, lab="Simulation", lw=0)
plot!(5dt:dt:5τ, t -> pdf(Exponential(τ), t), lw=2, lab="exp(-t/τ)/τ")
plot!(xlab="t", ylab="P(t)")
