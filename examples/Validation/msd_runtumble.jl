# # Mean squared displacement

#=
The mean squared displacement (MSD) of a population of run-tumble swimmers
with an average turn angle ``\theta``, should obey, in the absence of
rotational diffusion, the following equation (Taktikos *et al* (PLoS ONE 2012)):
```math
\text{MSD}(t) = 2v^2\tilde{\tau}^2
\left( \dfrac{t}{\tilde{\tau}} - 1 + \exp(-t/\tilde{\tau}) \right)
```
where ``v`` is the microbe velocity, and ``\tilde{\tau} = \tau/(1-\cos\theta)``
is the correlation-corrected effective run length.

To test this, we will simulate run-tumble motility with different distributions
of `polar` reorientation angles and evaluate their MSDs.
For simplicity, since only the average angle matters rather than the entire
distribution, we will only allow microbes to perform rotations of angle `θ`
(and `-θ` for symmetry).
=#

using MicrobeAgents
using Plots

dt = 0.05 # s
L = 500.0 # μm
extent = (L,L,L)
space = ContinuousSpace(extent)

θs = [π/6, π/4, π/3, π/2, π]
α = cos.(θs)

U = 30.0 # μm/s
τ = 1.0 # s
turn_rate = 1 / τ

## we initialize a separate model for each different θ
models = map(_ -> StandardABM(Microbe{3}, space, dt; container=Vector), θs)
nmicrobes = 100
for (i,θ) in enumerate(θs)
    motility = RunTumble(speed=[U], polar=[θ,-θ])
    foreach(_ -> add_agent!(models[i]; motility, turn_rate), 1:nmicrobes)
end

nsteps = round(Int, 100τ / dt)
adata = [position]
adfs = [run!(model, nsteps; adata)[1] for model in models]

foreach(adf -> Analysis.unfold!(adf, extent), adfs)
MSD = hcat([Analysis.emsd(adf, :position_unfold)[2:end] for adf in adfs]...)

t = (1:nsteps).*dt
β = cos.(θs)
T = @. τ / (1-β')
s = t ./ T
D = @. U^2*T/3
MSD_theoretical = @. 6D*T * (s - 1 + exp(-s))
plot(
    xlab = "Δt / τ",
    ylab = "MSD / (Uτ)²",
    legend = :bottomright, legendtitle = "1-cosθ",
    scale = :log10,
    yticks = exp10.(-2:2:2),
    xticks = exp10.(-2:2)
)
## show simulation values only at selected lags
lags = round.(Int, exp10.(range(0, 3, length=20))) |> unique
scatter!(t[lags]./τ, MSD[lags,:]./(U*τ)^2,
    m=:x, ms=6, msw=2, lab=false, lc=axes(β,1)'
)
plot!(t./τ, MSD_theoretical./(U*τ)^2,
    lw=2, lab=round.(1 .- β,digits=2)', lc=axes(β,1)'
)
