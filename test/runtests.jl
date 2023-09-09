using MicrobeAgents
using Test

@testset "MicrobeAgents.jl" begin
    include("utils.jl")
    include("distributions.jl")
    include("motility.jl")
    include("model-creation.jl")
    include("model-stepping.jl")
    include("neighborlists.jl")
    include("spheres.jl")
    include("analysis.jl")
end
