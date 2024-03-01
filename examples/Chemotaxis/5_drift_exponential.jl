# # Chemotactic drift in exponential ramp with noisy sensing

#=
We will peform simulations of chemotaxis in an exponential concentration
ramp using the `Brumley` model of chemotaxis and measure their drift velocity
along the gradient.
The concentration field has the form
```math
C(x) = C_0\exp(x/λ)
```
and its gradient is
```math
\nabla C(x) = \dfrac{C_0}{\lambda}\exp(x/\lambda)\hat{\mathbf{x}}.
```

In the `Brumley` model, the microbe estimates the local gradient in the form of
the convective derivative of the concentration field (``U\nabla C``).
The measurement is affected by noise, represented by an intrinsic
noise term ``\sigma = \sqrt{3C/(\pi a D_c T^3)}`` (Mora & Wingreen (Phys Rev Lett 2010),
where ``a`` is the microbe radius, ``T`` the chemotactic sensory timescale
and ``D_c`` the thermal diffusivity of the chemoattractant compound)
and a pre-factor ``\Pi`` which represents the "chemotactic precision".
Higher values of ``\Pi`` imply higher levels of noise in chemotaxis pathway.

Therefore, with each measurement the bacterium samples from a Normal distribution
with mean ``U\nabla C`` and standard deviation ``\Pi\sigma``.
The measurement, ``\mathcal{M}``, determines the evolution of an internal state
``S`` which obeys an excitation-relaxation dynamics:
```math
\dot{S} = \kappa M -\dfrac{S}{\tau_M}.
```
This internal state, in turn, determines the tumbling rate.

We will study how ``\Pi`` affects the drift velocity of bacteria.
=#

using MicrobeAgents

function concentration_field(pos, model)
    x = first(pos)
    C0 = model.C0
    λ = model.λ
    concentration_field(x, C0, λ)
end
concentration_field(x, C0, λ) = C0*exp(x/λ)

function concentration_gradient(pos, model)
    x = first(pos)
    C0 = model.C0
    λ = model.λ
    concentration_gradient(x, C0, λ)
end
concentration_gradient(x, C0, λ) = SVector{3}(i==1 ? C0/λ*exp(x/λ) : 0.0 for i in 1:3)

## wide rectangular channel confined in the z direction
Lx, Ly, Lz = 6000, 3000, 100
space = ContinuousSpace((Lx,Ly,Lz); periodic=false)
dt = 0.1

C0 = 10.0
λ = Lx/2
chemoattractant = GenericChemoattractant{3,Float64}(;
    concentration_field,
    concentration_gradient
)
properties = Dict(
    :chemoattractant => chemoattractant,
    :C0 => C0,
    :λ => λ,
)
model = StandardABM(Brumley{3}, space, dt; properties, container=Vector)

## add n bacteria for each value of Π
## all of them initialized at x = 0
Πs = [1.0, 2.0, 5.0, 10.0, 25.0]
n = 200
for Π in Πs
    for i in 1:n
        pos = SVector{3}(0.0, rand()*Ly, rand()*Lz)
        add_agent!(pos, model; chemotactic_precision=Π)
    end
end

## also store Π for grouping later
adata = [position, velocity, :chemotactic_precision]
nsteps = 2000
adf, = run!(model, nsteps; adata)

## evaluate drift velocity along x direction
target_direction = SVector(1.0, 0.0, 0.0)
Analysis.driftvelocity_direction!(adf, target_direction)

## we now want to first average the drift velocities in each group of Π values
## and then average each group over time to obtain a mean drift
## the calculation is much easier to perform with the DataFrames package
using DataFrames, StatsBase
gdf = groupby(adf, :chemotactic_precision)
drift_velocities = zeros(1+nsteps, length(Πs))
for (i,g) in enumerate(gdf)
    for h in groupby(g, :id)
        drift_velocities[:,i] .+= h.drift_direction
    end
    drift_velocities[:,i] ./= n
end
avg_drift = vec(mean(drift_velocities; dims=1))

## we apply a smoothing function to remove some noise for better visualization
function moving_average(y::AbstractVector, m::Integer)
    @assert isodd(m)
    out = similar(y)
    R = CartesianIndices(y)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += y[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end
smoothing_window = 51
smooth_drift_velocities = mapslices(
    y -> moving_average(y, smoothing_window), drift_velocities;
    dims = 1
)

using Plots
t = axes(smooth_drift_velocities, 1) .* dt
plot(layout=(1,2), size=(760,440),
    plot(t, smooth_drift_velocities, lab=Πs',
        xlab="t (s)", ylab="drift (μm/s)"
    ),
    plot(Πs, avg_drift, lw=4, m=:c, ms=8, lab=false,
        xlab="Π", ylab="avg drift (μm/s)"
    )
)

#=
Increasing Π significantly reduces the ability of bacteria to drift
along the gradient.
=#
