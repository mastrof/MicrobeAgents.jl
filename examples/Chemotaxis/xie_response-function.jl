using MicrobeAgents
using Plots

θ(a,b) = a>b ? 1.0 : 0.0
function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    t₁ = model.t₁
    t₂ = model.t₂
    dt = model.timestep
    t = model.t * dt
    concentration_field(t, C₀, C₁, t₁, t₂)
end
concentration_field(t,C₀,C₁,t₁,t₂) = C₀+C₁*θ(t,t₁)*(1-θ(t,t₂))

timestep = 0.1 # s
space = ContinuousSpace(ntuple(_ -> 500.0, 3)) # μm
C₀ = 0.01 # μM
C₁ = 5.0-C₀ # μM
T = 60.0 # s
t₁ = 20.0 # s
t₂ = 40.0 # s
properties = Dict(
    :concentration_field => concentration_field,
    :C₀ => C₀,
    :C₁ => C₁,
    :t₁ => t₁,
    :t₂ => t₂,
)

model_step!(model) = model.t += 1
model = UnremovableABM(Xie{3}, space, timestep; properties, model_step!)
add_agent!(model; turn_rate_forward=0, motility=RunReverseFlick(motile_state=MotileState(Forward)))
add_agent!(model; turn_rate_backward=0, motility=RunReverseFlick(motile_state=MotileState(Backward)))

nsteps = round(Int, T/timestep)
β(a) = a.motility.state == Forward ? a.gain_forward : a.gain_backward
state(a::Xie) = max(1 + β(a)*a.state, 0)
adata = [state, :state_m, :state_z]
adf, = run!(model, nsteps; adata)

S = vectorize_adf_measurement(adf, :state)
m = (vectorize_adf_measurement(adf, :state_m))[:,1] # take only fw
z = (vectorize_adf_measurement(adf, :state_z))[:,1] # take only fw

# response vs time for fw and bw modes
begin
    _green = palette(:default)[3]
    plot()
    x = (0:timestep:T) .- t₁
    plot!(
        x, S,
        lw=1.5, lab=["Forward" "Backward"]
    )
    plot!(ylims=(-0.1,4.5), ylab="Response", xlab="time (s)")
    plot!(twinx(),
        x, t -> concentration_field(t.+t₁,C₀,C₁,t₁,t₂),
        ls=:dash, lw=1.5, lc=_green, lab=false,
        tickfontcolor=_green,
        ylab="C (μM)", guidefontcolor=_green
    )
end

# methylation and dephosphorylation
begin
    x = (0:timestep:T) .- t₁
    τ_m = model[1].adaptation_time_m
    τ_z = model[1].adaptation_time_z
    M = m ./ τ_m
    Z = z ./ τ_z
    R = M .- Z
    plot(
        x, [M Z R],
        lw=2,
        lab=["m/τ_m" "z/τ_z" "m/τ_m - z/τ_z"],
        xlab="time (s)"
    )
end
