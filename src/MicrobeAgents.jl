module MicrobeAgents

using Agents
using CellListMap.PeriodicSystems
using Distributions
using LinearAlgebra
using Random
using Rotations
using StaticArrays

using DSP, AxisArrays, OffsetArrays # for acf and msd


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
- `pos::NTuple{D,Float64}` position of the microbe
- `vel::NTuple{D,Float64}` velocity of the microbe
- `motility::AbstractMotility` motile pattern of the microbe
- `turn_rate::Real` average reorientation rate of the microbe
- `rotational_diffusivity::Real` coefficient of brownian rotational diffusion
- `radius::Real` equivalent spherical radius of the microbe
- `state::Real` generic variable for a scalar internal state
"""
abstract type AbstractMicrobe{D} <: AbstractAgent where D end

include("utils.jl")
include("distributions.jl")
include("motility.jl")
include("rotations.jl")
include("neighborlists.jl")

include("microbes.jl")
include("microbe_step.jl")
include("model.jl")

# implementations of chemotactic models
include("chemotaxis/brown-berg.jl")
include("chemotaxis/brumley.jl")
include("chemotaxis/celani.jl")
include("chemotaxis/xie.jl")

using Agents.Pathfinding
include("pathfinder.jl")

# analysis routines
include("analysis/postprocess.jl")
include("analysis/msd.jl")
include("analysis/correlation_functions.jl")
include("analysis/drift.jl")
include("analysis/turn_detection.jl")

# submodules
include("submodules/SphericalSurfaces.jl")

end
