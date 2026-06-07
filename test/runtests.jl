using MicrobeAgents
using Test

@testset "MicrobeAgents.jl" begin
    include("utils.jl")
    include("motility.jl")
    include("model-creation.jl")
    include("model-stepping.jl")
    include("analysis.jl")
end
