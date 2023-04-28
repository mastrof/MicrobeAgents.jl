push!(LOAD_PATH, "../src/")
using Documenter, MicrobeAgents

makedocs(
    sitename = "MicrobeAgents.jl",
    modules = [MicrobeAgents],
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Validation" => "validation.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(;
    repo = "github.com/mastrof/MicrobeAgents.jl"
)