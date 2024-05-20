cd(@__DIR__)
ENV["JULIA_DEBUG"] = "Documenter"
using MicrobeAgents
using Documenter
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
import Literate
using Plots

# use Literate to compile examples
indir_base = joinpath(@__DIR__, "..", "examples")
sections = ("RandomWalks", "Chemotaxis")
outdir_base = joinpath(@__DIR__, "src", "examples")
indir = Dict(
    section => joinpath(indir_base, section)
    for section in sections
)
outdir = Dict(
    section => joinpath(outdir_base, section)
    for section in sections
)
rm(outdir_base; force=true, recursive=true) # clean up previous examples
mkpath(outdir_base)
for section in sections
    mkpath(outdir[section])
    for file in readdir(indir[section])
        Literate.markdown(joinpath(indir[section], file), outdir[section]; credit=false)
    end
end

# use Literate to compile validation scripts
valdir_in = joinpath(@__DIR__, "..", "examples/Validation")
valdir_out = joinpath(@__DIR__, "src", "validation")
rm(valdir_out; force=true, recursive=true)
mkpath(valdir_out)
for file in readdir(valdir_in)
    Literate.markdown(joinpath(valdir_in, file), valdir_out; credit=false)
end


# convert camelcase directory names to space-separated section names
function namify(s)
    indices = findall(isuppercase, s)
    if length(indices) <= 1
        return s
    else
        s1 = s[indices[1] : indices[2]-1]
        s2 = lowercase.(s[indices[2] : end])
        return join((s1, namify(s2)), " ")
    end
end

pages = [
    "Home" => "index.md",
    "Introduction" => "introduction.md",
    "Examples" => [
        [namify(section) => [joinpath.("examples", section, readdir(outdir[section]))...]
         for section in sections]...
    ],
    "Validation" => joinpath.("validation", readdir(valdir_out)),
    "API" => "api.md"
]

makedocs(
    sitename = "MicrobeAgents.jl",
    authors = "Riccardo Foffi",
    modules = [MicrobeAgents],
    pages = pages,
    expandfirst = ["index.md"],
    format = Documenter.HTML(
        prettyurls = CI,
    ),
    warnonly = true,
)

if CI
    deploydocs(;
        repo = "github.com/mastrof/MicrobeAgents.jl.git",
        push_preview = true
    )
end
