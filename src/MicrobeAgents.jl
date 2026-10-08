module MicrobeAgents

using Agents
export abmproperties, abmrng, abmspace, abmscheduler, abmtime, spacesize
export ABM, StandardABM, ContinuousSpace
export add_agent!, add_agent_own_pos!
export move_agent!, walk!, run!

using Distributions
using LightSumTypes
using LinearAlgebra
using Random
using Quaternions
using StaticArrays
using StatsBase
export SVector


export AbstractMicrobe
"""
    AbstractMicrobe{D,N} <: AbstractAgent where {D<:Integer,N}
All microbe types in MicrobeAgents.jl simulations must be instances
of user-defined types that are subtypes of `AbstractMicrobe`.
    YourMicrobeType{D,N} <: AbstractMicrobe{D,N}
The parameter `D` defines the dimensionality of the space in which the
microbe type lives (1, 2 and 3 are supported); `N` is the number of motile
states of its motility.

All microbe types *must* have at least the following fields:
- `id::Int` id of the microbe (used internally by Agents.jl)
- `pos::SVector{D,Float64}` position of the microbe
- `vel::SVector{D,Float64}` velocity of the microbe
- `speed::Real` base speed of the microbe (see `speed`)
- `motility::Motility{N}` motile pattern of the microbe
- `rotational_diffusivity::Real` coefficient of brownian rotational diffusion
- `radius::Real` equivalent spherical radius of the microbe

Optionally, define a method `behaviors(m)` returning a `Tuple` or `NamedTuple`
of behaviors (default: `()`).
"""
abstract type AbstractMicrobe{D,N} <: AbstractAgent where {D,N} end

include("api.jl")
include("utils.jl")
include("spherical_distribution.jl")
include("motility.jl")
include("rotations.jl")
include("fields.jl")
include("behaviors.jl")

include("microbes.jl")
include("microbe_step.jl")
include("model.jl")

# implementations of chemotactic models as behaviors
include("chemotaxis/sensing.jl")
include("chemotaxis/brown-berg.jl")
include("chemotaxis/brumley.jl")
include("chemotaxis/celani.jl")
include("chemotaxis/xie.jl")
include("chemotaxis/son-menolascina.jl")

# pathfinding
using Agents.Pathfinding
include("pathfinder.jl")

# submodules
include("submodules/Analysis/Analysis.jl")

end
