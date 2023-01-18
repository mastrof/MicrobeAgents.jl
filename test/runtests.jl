using MicrobeAgents
using Test

@testset "MicrobeAgents.jl" begin
    include("utils.jl")
    include("distributions.jl")
    include("motility.jl")
    include("microbe-creation.jl")
    include("model.jl")
    include("analysis.jl")
end
