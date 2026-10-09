export AbstractChemicalField, ChemicalField, chemicalfield
export concentration, gradient, time_derivative, diffusivity

"""
    AbstractChemicalField{D}
Abstract type for chemical fields.
Requires dimensionality (`D`) to be specified.
Number type is always assumed to be `Float64`.

The interface is defined by four functions:
- `concentration`: returns the function for the concentration field
- `gradient`: returns the function for the concentration gradient
- `time_derivative`: returns the function for the concentration ramp
- `diffusivity`: returns the diffusivity of the chemical
"""
abstract type AbstractChemicalField{D} end

# per-step memoization of field quantities for the microbe being stepped;
# `id == 0` means no microbe is being stepped.
# Each field has a slot whose `valid` bitmask records the quantities computed so
# far, so a reset only clears the masks and stale values are never read.
mutable struct FieldSlot{D}
    valid::UInt8
    concentration::Float64
    gradient::SVector{D,Float64}
    time_derivative::Float64
    diffusivity::Float64
end
FieldSlot{D}() where {D} = FieldSlot{D}(0x00, 0.0, zero(SVector{D,Float64}), 0.0, 0.0)

# `keys[i]`, `fields[i]` and `slots[i]` describe the i-th field (`:chemicalfield` first);
# fields are looked up by name with a short linear scan of `keys`
mutable struct FieldCache{D}
    id::Int
    keys::Vector{Symbol}
    fields::Vector{AbstractChemicalField{D}}
    slots::Vector{FieldSlot{D}}
end

# collect every `AbstractChemicalField` property, default field first
function FieldCache{D}(props) where {D}
    default = getproperty(props, :chemicalfield)
    default isa AbstractChemicalField{D} || throw(ArgumentError(
        "property `:chemicalfield` must be an `AbstractChemicalField{$D}`, got $(typeof(default))"))
    for k in keys(props)
        v = getproperty(props, k)
        v isa AbstractChemicalField && !(v isa AbstractChemicalField{D}) && throw(ArgumentError(
            "chemical field `:$k` has the wrong dimension: expected `AbstractChemicalField{$D}`, got $(typeof(v))"))
    end
    ks = [k for k in keys(props) if k !== :chemicalfield && getproperty(props, k) isa AbstractChemicalField]
    pushfirst!(ks, :chemicalfield)
    fields = AbstractChemicalField{D}[getproperty(props, k) for k in ks]
    FieldCache{D}(0, ks, fields, [FieldSlot{D}() for _ in ks])
end

const _CONCENTRATION_BIT = 0x01
const _GRADIENT_BIT = 0x02
const _TIME_DERIVATIVE_BIT = 0x04
const _DIFFUSIVITY_BIT = 0x08

field_cache(model::ABM) = abmproperties(model).field_cache

function reset_field_cache!(model::ABM, microbe::AbstractMicrobe)
    c = field_cache(model)
    c.id = microbe.id
    for s in c.slots
        s.valid = 0x00
    end
    return nothing
end
"""
    MicrobeAgents.invalidate_field_cache!(model)
Discard cached field values for the current step. Call it after a behavior
changes the position of the microbe, so later behaviors see fresh values.
"""
invalidate_field_cache!(model::ABM) = (field_cache(model).id = 0; nothing)

@noinline function _unknown_field(cache::FieldCache, key::Symbol)
    throw(ArgumentError(
        "model has no chemical field `:$key`; available fields: " *
        join((":$k" for k in cache.keys), ", ")))
end

# position of `key` in the (short) vector of field names
@inline function _index(cache::FieldCache, key::Symbol)
    ks = cache.keys
    @inbounds for i in eachindex(ks)
        ks[i] === key && return i
    end
    return _unknown_field(cache, key)
end

"""
    MicrobeAgents.check_field(model, key::Symbol)
Throw an `ArgumentError` unless `model` has a chemical field stored under the
property `key`.
"""
check_field(model::ABM, key::Symbol) = (_index(field_cache(model), key); nothing)

# `compute(field, microbe, model)` evaluates the quantity (a non-capturing function,
# so that no closure is allocated when it is passed to a non-inlined method).
# The default field is read from the concretely typed model property, so that
# `compute` dispatches statically on the single-field path (slot 1 is always the default).
@inline function _cached(compute::F, name::Symbol, bit::UInt8, microbe, model, key) where {F}
    cache = field_cache(model)
    if key === :chemicalfield
        return _cached_slot(compute, chemicalfield(model), cache, 1, name, bit, microbe, model)
    end
    i = _index(cache, key)
    return @inline _cached_slot(compute, cache.fields[i], cache, i, name, bit, microbe, model)
end

@inline function _cached_slot(compute::F, field, cache, i, name, bit, microbe, model) where {F}
    cache.id == microbe.id || return compute(field, microbe, model)
    slot = @inbounds cache.slots[i]
    slot.valid & bit == bit && return getfield(slot, name)
    v = compute(field, microbe, model)
    setfield!(slot, name, v)
    slot.valid |= bit
    return v
end

"""
    concentration(microbe, model[, key])
Concentration at the position of `microbe` in the chemical field `key`
(default `:chemicalfield`). Memoized per step.
"""
function concentration(microbe::AbstractMicrobe{D,N}, model::ABM, key::Symbol) where {D,N}
    _cached(:concentration, _CONCENTRATION_BIT, microbe, model, key) do f, microbe, model
        concentration(f)(microbe, model)::Float64
    end
end
"""
    gradient(microbe, model[, key])
Concentration gradient for `microbe` in the chemical field `key` (default `:chemicalfield`).
"""
function gradient(microbe::AbstractMicrobe{D,N}, model::ABM, key::Symbol) where {D,N}
    _cached(:gradient, _GRADIENT_BIT, microbe, model, key) do f, microbe, model
        gradient(f)(microbe, model)::SVector{D,Float64}
    end
end
"""
    time_derivative(microbe, model[, key])
Time derivative of the concentration for `microbe` in the chemical field `key` (default `:chemicalfield`).
"""
function time_derivative(microbe::AbstractMicrobe{D,N}, model::ABM, key::Symbol) where {D,N}
    _cached(:time_derivative, _TIME_DERIVATIVE_BIT, microbe, model, key) do f, microbe, model
        time_derivative(f)(microbe, model)::Float64
    end
end
"""
    diffusivity(microbe, model[, key])
Diffusivity of the chemical `key` (default `:chemicalfield`) at the position of `microbe`.
"""
function diffusivity(microbe::AbstractMicrobe{D,N}, model::ABM, key::Symbol) where {D,N}
    _cached(:diffusivity, _DIFFUSIVITY_BIT, microbe, model, key) do f, microbe, model
        diffusivity(f)(microbe, model)::Float64
    end
end
concentration(microbe::AbstractMicrobe, model::ABM) = concentration(microbe, model, :chemicalfield)
gradient(microbe::AbstractMicrobe, model::ABM) = gradient(microbe, model, :chemicalfield)
time_derivative(microbe::AbstractMicrobe, model::ABM) = time_derivative(microbe, model, :chemicalfield)
diffusivity(microbe::AbstractMicrobe, model::ABM) = diffusivity(microbe, model, :chemicalfield)

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
Returns the thermal diffusivity of the default chemical field.
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
