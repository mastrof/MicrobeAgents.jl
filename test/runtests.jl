using MicrobeAgents
using Test

@testset "MicrobeAgents.jl" begin
    include("utils.jl")
    include("behaviors.jl")
    include("field-cache.jl")
    include("motility.jl")
    include("model-creation.jl")
    include("model-stepping.jl")
    # include("chemotaxis.jl") # re-enabled as behaviors in later tasks
    include("analysis.jl")
end
