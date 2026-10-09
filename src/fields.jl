export AbstractChemicalField, ChemicalField, chemicalfield
export concentration, gradient, time_derivative, diffusivity

"""
    AbstractChemicalField{D}
Abstract type for chemoattractants.
Requires dimensionality (`D`) to be specified.
Number type is always assumed to be `Float64`.

The interface is defined by five core functions:
- `chemicalfield`: returns the default chemical field object
- `concentration`: returns the function for the concentration field
- `gradient`: returns the function for the concentration gradient
- `time_derivative`: returns the function for the concentration ramp
- `diffusivity`: returns the thermal diffusivity of the chemoattractant
"""
abstract type AbstractChemicalField{D} end

# per-step memoization of field quantities for the microbe being stepped;
# `id == 0` means no microbe is being stepped.
# `valid` is a bitmask of the quantities computed so far, so a reset only
# clears the mask and the stale values are never read.
mutable struct FieldCache{D}
    id::Int
    valid::UInt8
    concentration::Float64
    gradient::SVector{D,Float64}
    time_derivative::Float64
    diffusivity::Float64
end
FieldCache{D}() where {D} = FieldCache{D}(0, 0x00, 0.0, zero(SVector{D,Float64}), 0.0, 0.0)

const _CONCENTRATION_BIT = 0x01
const _GRADIENT_BIT = 0x02
const _TIME_DERIVATIVE_BIT = 0x04
const _DIFFUSIVITY_BIT = 0x08

field_cache(model::ABM) = abmproperties(model).field_cache

function reset_field_cache!(model::ABM, microbe::AbstractMicrobe)
    c = field_cache(model)
    c.id = microbe.id
    c.valid = 0x00
    return nothing
end
"""
    MicrobeAgents.invalidate_field_cache!(model)
Discard cached field values for the current step. Call it after a behavior
changes the position of the microbe, so later behaviors see fresh values.
"""
invalidate_field_cache!(model::ABM) = (field_cache(model).id = 0; nothing)

@inline function _cached(compute::F, name::Symbol, bit::UInt8, microbe, model) where {F}
    cache = field_cache(model)
    cache.id == microbe.id || return compute()
    cache.valid & bit == bit && return getfield(cache, name)
    v = compute()
    setfield!(cache, name, v)
    cache.valid |= bit
    return v
end

function concentration(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:concentration, _CONCENTRATION_BIT, microbe, model) do
        concentration(chemicalfield(model))(microbe, model)::Float64
    end
end
function gradient(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:gradient, _GRADIENT_BIT, microbe, model) do
        gradient(chemicalfield(model))(microbe, model)::SVector{D,Float64}
    end
end
function time_derivative(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:time_derivative, _TIME_DERIVATIVE_BIT, microbe, model) do
        time_derivative(chemicalfield(model))(microbe, model)::Float64
    end
end
function diffusivity(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:diffusivity, _DIFFUSIVITY_BIT, microbe, model) do
        diffusivity(chemicalfield(model))(microbe, model)::Float64
    end
end

"""
    chemicalfield(model)
Returns the default chemical field of `model`.
"""
chemicalfield(model::ABM) = model.chemicalfield
"""
    concentration(model)
Returns the function `f` that defines the concentration field.
The returned function has signature `f(pos, model)` and returns a scalar.
"""
concentration(model::ABM) = concentration(chemicalfield(model))
"""
    gradient(model)
Returns the function `f` that defines the gradient of the concentration field.
The returned function has signature `f(pos, model)` and returns a `SVector`
with the same dimensionality as the microbe position `pos`.
"""
gradient(model::ABM) = gradient(chemicalfield(model))
"""
    time_derivative(model)
Returns the function `f` that defines the time derivative of the concentration field.
The returned function has signature `f(pos, model)` and returns a scalar.
"""
time_derivative(model::ABM) = time_derivative(chemicalfield(model))
"""
    diffusivity(model)
Returns the thermal diffusivity of the chemoattractant compound.
"""
diffusivity(model::ABM) = diffusivity(chemicalfield(model))
concentration(c::AbstractChemicalField) = c.concentration_field
gradient(c::AbstractChemicalField) = c.concentration_gradient
time_derivative(c::AbstractChemicalField) = c.concentration_ramp
diffusivity(c::AbstractChemicalField) = c.diffusivity

"""
    ChemicalField{D} <: AbstractChemicalField{D}
Type for a generic chemical field.
Field, gradient and ramp default to 0 everywhere in the domain.
Diffusivity defaults to 608 μm²/s everywhere in the domain. (Do not use 0 here, it may mess up some calculations)
"""
struct ChemicalField{D} <: AbstractChemicalField{D}
    concentration_field::Function
    concentration_gradient::Function
    concentration_ramp::Function
    diffusivity::Function
end
function ChemicalField{D}(;
    concentration_field = (::AbstractMicrobe, ::ABM) -> zero(Float64), # μM
    concentration_gradient = (::AbstractMicrobe, ::ABM) -> zero(SVector{D,Float64}), # μM/μm
    concentration_ramp = (::AbstractMicrobe, ::ABM) -> zero(Float64), # μM/s
    diffusivity = (::AbstractMicrobe, ::ABM) -> Float64(608), # μm²/s
) where D
    ChemicalField{D}(
        concentration_field, concentration_gradient, concentration_ramp, diffusivity
    )
end
