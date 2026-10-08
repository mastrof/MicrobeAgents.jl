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
    AbstractMicrobe{D} <: AbstractAgent where {D<:Integer}
All microbe types in MicrobeAgents.jl simulations must be instances
of user-defined types that are subtypes of `AbstractMicrobe`.
    YourMicrobeType{D} <: AbstractMicrobe{D}
The parameter `D` defines the dimensionality of the space in which the
microbe type lives (1, 2 and 3 are supported).

All microbe types *must* have at least the following fields:
- `id::Int` id of the microbe (used internally by Agents.jl)
- `pos::SVectpr{D,Float64}` position of the microbe
- `vel::SVector{D,Float64}` velocity of the microbe
- `speed::Real` speed of the microbe
- `motility::AbstractMotility` motile pattern of the microbe
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

# implementations of chemotactic models
"""
Conversion factor (1/√(number of molecules) --> 1/√(moles)) used
in the evaluation of chemotactic sensing noise.
"""
global const CONV_NOISE::Float64 = 0.04075
# include("chemotaxis/brown-berg.jl") # re-enabled as behaviors in later tasks
# include("chemotaxis/brumley.jl") # re-enabled as behaviors in later tasks
# include("chemotaxis/celani.jl") # re-enabled as behaviors in later tasks
# include("chemotaxis/xie.jl") # re-enabled as behaviors in later tasks
# include("chemotaxis/son-menolascina.jl") # re-enabled as behaviors in later tasks

# pathfinding
using Agents.Pathfinding
include("pathfinder.jl")

# submodules
include("submodules/Analysis/Analysis.jl")

end
