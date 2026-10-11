# # Rectification of motile bacteria by a wall of funnels

#=
A wall of funnel-shaped obstacles breaks the left-right symmetry of the motion
of swimming bacteria: those coming from the left are guided through the openings
by the slanted arms, while those coming from the right are trapped in the
concave vertices.
The result is a net accumulation of bacteria on the right side of the wall
(Galajda et al., *J. Bacteriol.* **189**, 8704 (2007)).

Each funnel is a pair of capsules joined at their left end, forming a `<` that
opens toward the right. Bacteria slide along the surfaces and are reflected by
the lateral walls of the chamber, while the vertical direction is periodic.
=#

using MicrobeAgents
using Random
using Plots

dt = 0.1 # s
L = 300.0 # μm
space = ContinuousSpace((L, L); periodic = (false, true))

# ## funnel wall
nf = 6 # number of funnels
p = L / nf # funnel period
len = 30.0 # arm length
α = deg2rad(40) # arm angle with respect to the horizontal
r = 1.0 # arm half-thickness
## the gap between the tips of adjacent funnels is p - 2len*sin(α) - 2r

function funnel(vertex)
    tip(s) = vertex + len * SVector(cos(α), s * sin(α))
    [ObstacleCapsule(vertex, tip(+1), r), ObstacleCapsule(vertex, tip(-1), r)]
end
x0 = L / 2 - len * cos(α) / 2
arms = reduce(vcat, [funnel(SVector(x0, (k + 1 // 2) * p)) for k in 0:(nf-1)])

obstacles = (
    Obstacles(space, arms; contact = Slide()),
    Walls(space; contact = Reflect()),
)
model = StandardABM(Microbe{2}, space, dt; rng = Xoshiro(1), container = Vector,
    properties = Dict(:obstacles => obstacles))

# ## bacteria, initially uniformly distributed on both sides
motility = RunTumble(run_speed = [20.0], run_duration = 1.0, angle = Isotropic(2))
for _ in 1:300
    pos = random_free_position(model; radius = 0.5)
    add_agent!(pos, model; motility, radius = 0.5, rotational_diffusivity = 0.035)
end

nsteps = 15_000
adf, _ = run!(model, nsteps; adata = [position], when = 0:100:nsteps)
X = Analysis.adf_to_matrix(adf, :position)
xs = first.(X) # rows: times, columns: bacteria
right = vec(sum(xs .> L / 2; dims = 2)) ./ size(xs, 2)

# ## results
function capsule(c)
    pts = vcat(
        [c.a + c.radius * SVector(cos(t), sin(t)) for t in range(π / 2, 3π / 2, 16)],
        [c.b + c.radius * SVector(cos(t), sin(t)) for t in range(-π / 2, π / 2, 16)],
    )
    Shape(first.(pts), last.(pts))
end
## fraction of bacteria on the right of the wall as a function of time
t = (0:100:nsteps) .* dt
p1 = plot(t, right, xlabel = "time (s)", ylabel = "fraction on the right",
    ylims = (0, 1), leg = false, lw = 2)
hline!(p1, [0.5], c = :gray, ls = :dash)
## final configuration
p2 = plot(xlims = (0, L), ylims = (0, L), ratio = 1, axis = false, grid = false, leg = false)
foreach(c -> plot!(p2, capsule(c), c = :gray, lw = 0), arms)
scatter!(p2, xs[end, :], last.(X)[end, :], ms = 2, c = :green, msw = 0)
plot(p1, p2, size = (900, 400))
