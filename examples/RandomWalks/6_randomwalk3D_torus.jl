# # Random walk around a custom 3D obstacle

#=
Obstacles of arbitrary shape can be defined from their signed distance
function with `ObstacleSDF`. Here we place a torus at the center of a
three-dimensional box with reflecting walls; bacteria slide along its surface.

The function must return a negative value inside the solid and must not
overestimate the distance to the surface (the torus distance is exact).
The center and radius of a ball enclosing the solid let the obstacle be
indexed spatially. The torus is not convex, so we leave `convex = false`
(the default): collisions are then resolved by sphere tracing.
=#

using MicrobeAgents
using Random
using Plots

dt = 0.02 # s
L = 200.0 # μm
space = ContinuousSpace((L, L, L); periodic = false)

## torus of major radius R and minor radius r, with its axis along z
R, r = 50.0, 20.0
c = SVector(L / 2, L / 2, L / 2)
function torus(x)
    p = x - c
    hypot(hypot(p[1], p[2]) - R, p[3]) - r
end
obstacle = ObstacleSDF(torus, c, R + r)

obstacles = (
    Obstacles(space, [obstacle]; contact = Slide(; align=true)),
    Walls(space; contact=Halt()),
)
rng = Xoshiro(1)
model = StandardABM(Microbe{3}, space, dt; rng, container = Vector,
    properties = Dict(:obstacles => obstacles))

motility = RunTumble(run_speed = [20.0], run_duration = Inf, angle = Isotropic(3))
for _ in 1:1
    pos = SVector{3}(L/2-R-r-20, L/2, L/2)
    vel = SVector{3}(1.0, 0.01, 0.1)
    vel = vel / sqrt(sum(abs2.(vel)))
    add_agent!(pos, model; vel, motility, radius = 0.5, rotational_diffusivity = 0)
end
adf, _ = run!(model, 3000; adata = [position], when = 2)
traj = Analysis.adf_to_matrix(adf, :position)
x, y, z = first.(traj), getindex.(traj, 2), last.(traj)
## bacteria never enter the torus (distance from its surface ≥ their radius)
# @assert minimum(torus.(traj)) ≥ 0.5 - 1e-6

# ## results
θ = range(0, 2π, 60)
X = [c[1] + (R + r * cos(v)) * cos(u) for v in θ, u in θ]
Y = [c[2] + (R + r * cos(v)) * sin(u) for v in θ, u in θ]
Z = [c[3] + r * sin(v) for v in θ, u in θ]
surface(X, Y, Z; c = :gray, alpha = 0.3, cbar = false, leg = false, ratio = 1,
    xlims = (0, L), ylims = (0, L), zlims = (0, L), xlab = "x", ylab = "y", zlab = "z")
plot!(x, y, z; lw = 1)
