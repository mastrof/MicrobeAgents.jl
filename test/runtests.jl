using MicrobeAgents
using Test

@testset "MicrobeAgents.jl" begin
    include("utils.jl")
    include("behaviors.jl")
    include("fields.jl")
    include("regression.jl")
    include("motility.jl")
    include("model-creation.jl")
    include("model-stepping.jl")
    include("custom-types.jl")
    include("obstacles.jl")
    include("chemotaxis.jl")
    include("analysis.jl")
    include("quality.jl")
end
