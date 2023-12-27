cd(@__DIR__)
using MicrobeAgents
import Literate

pages = [
    "Home" => "index.md",
    #"Tutorial" => ["firststeps.md", "randomwalks.md", "chemotaxis.md"],
    #"Validation" => "validation.md",
    #"API" => "api.md"
]

indir = joinpath(@__DIR__, "..", "examples", "RandomWalks")
outdir = joinpath(@__DIR__, "src", "examples")
rm(outdir; force=true, recursive=true) # clean up previous examples
mkpath(outdir)
toskip = ()
for file in readdir(indir)
    file in toskip && continue
    Literate.markdown(joinpath(indir, file), outdir; credit=false)
end

using Documenter
ENV["JULIA_DEBUG"] = "Documenter"
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

makedocs(
    sitename = "MicrobeAgents.jl",
    authors = "Riccardo Foffi",
    modules = [MicrobeAgents],
    pages = pages,
    format = Documenter.HTML(
        prettyurls = CI,
    ),
    warnonly = true,
)

if CI
    deploydocs(;
        repo = "github.com/mastrof/MicrobeAgents.jl.git",
        target = "build",
        push_preview = true
    )
end
