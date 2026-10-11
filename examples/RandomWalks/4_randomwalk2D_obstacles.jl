# # Random walk among obstacles

#=
Solid obstacles are described by signed distance functions and grouped in
`Obstacles` sets; each set defines what happens when a microbe touches its
surface (its contact rule). Sets are passed as a `Tuple` in the `:obstacles`
model property, and the domain boundaries can be made solid with `Walls`.

Here we generate a random packing of disks, on which bacteria slide,
inside a box with reflecting walls.
Unlike the pathfinder approach, surfaces are exact (no discretization).
=#

using MicrobeAgents
using Random
using Plots

dt = 0.1 # s
L = 300.0 # μm
space = ContinuousSpace((L, L); periodic = false)

## random sequential addition of non-overlapping disks
function disk_packing(rng, L, φ, rmin, rmax; maxtries = 100_000)
    disks = ObstacleSphere{2}[]
    area = 0.0
    for _ in 1:maxtries
        area ≥ φ * L^2 && break
        r = rmin + (rmax - rmin) * rand(rng)
        c = SVector(rand(rng) * L, rand(rng) * L)
        any(d -> sdf(d, c) < r, disks) && continue
        push!(disks, ObstacleSphere(c, r))
        area += π * r^2
    end
    disks
end
rng = Xoshiro(1)
disks = disk_packing(rng, L, 0.4, 10.0, 20.0)

obstacles = (
    Obstacles(space, disks; contact = Slide()),
    Walls(space; contact = Reflect()),
)
model = StandardABM(Microbe{2}, space, dt; rng, container = Vector,
    properties = Dict(:obstacles => obstacles))

motility = RunTumble(run_speed = [20.0], run_duration = 1.0, angle = Isotropic(2))
for _ in 1:10
    ## sample initial positions outside the obstacles
    pos = random_free_position(model; radius = 0.5)
    add_agent!(pos, model; motility, radius = 0.5, rotational_diffusivity = 0.035)
end
adf, _ = run!(model, 1000; adata = [position])
traj = Analysis.adf_to_matrix(adf, :position)

circle(d) = Shape(
    [d.center[1] + d.radius * cos(t) for t in range(0, 2π, 64)],
    [d.center[2] + d.radius * sin(t) for t in range(0, 2π, 64)],
)
plot(xlims = (0, L), ylims = (0, L), ratio = 1, axis = false, grid = false, leg = false)
foreach(d -> plot!(circle(d), c = :gray, lw = 0), disks)
plot!(first.(traj), last.(traj))
