using MicrobeAgents
using Plots
using Random

function concentration_field(pos, model)
    C = model.C
    σ = model.σ
    p₀ = model.p₀
    concentration_field(pos, p₀, C, σ)
end
concentration_field(pos, p₀, C, σ) = C * exp(-sum(abs2.(pos.-p₀))/(2*σ^2))

timestep = 0.1 # s
extent = ntuple(_ -> 1000.0, 2) # μm
p₀ = extent./2 # μm
C = 500.0 # μM
σ = 25.0 # μm
properties = Dict(
    :concentration_field => concentration_field,
    :C => C,
    :σ => σ,
    :p₀ => p₀,
)

rng = MersenneTwister(12)
model = ABM(Xie{2}, extent, timestep; rng, properties, periodic=false)
foreach(_ -> add_agent!(model; chemotactic_precision=6.0), 1:300)

nsteps = 5000
adata = [:pos]
adf, = run!(model, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
plot(
    first.(traj)[end-100:end,:],
    last.(traj)[end-100:end,:],
    lab=false, lims=(0,1000)
)