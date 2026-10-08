export AbstractChemoattractant, GenericChemoattractant, chemoattractant
export concentration, gradient, time_derivative, chemoattractant_diffusivity

"""
    AbstractChemoattractant{D}
Abstract type for chemoattractants.
Requires dimensionality (`D`) to be specified.
Number type is always assumed to be `Float64`.

The interface is defined by five core functions:
- `chemoattractant`: returns the chemoattractant object
- `concentration`: returns the function for the concentration field
- `gradient`: returns the function for the concentration gradient
- `time_derivative`: returns the function for the concentration ramp
- `chemoattractant_diffusivity`: returns the thermal diffusivity of the chemoattractant
"""
abstract type AbstractChemoattractant{D} end

# per-step memoization of field quantities for the microbe being stepped;
# `id == 0` means no microbe is being stepped
mutable struct FieldCache{D}
    id::Int
    concentration::Union{Nothing,Float64}
    gradient::Union{Nothing,SVector{D,Float64}}
    time_derivative::Union{Nothing,Float64}
    chemoattractant_diffusivity::Union{Nothing,Float64}
end
FieldCache{D}() where {D} = FieldCache{D}(0, nothing, nothing, nothing, nothing)

field_cache(model::ABM) = abmproperties(model).field_cache

function reset_field_cache!(model::ABM, microbe::AbstractMicrobe)
    c = field_cache(model)
    c.id = microbe.id
    c.concentration = nothing
    c.gradient = nothing
    c.time_derivative = nothing
    c.chemoattractant_diffusivity = nothing
    return nothing
end
"""
    MicrobeAgents.invalidate_field_cache!(model)
Discard cached field values for the current step. Call it after a behavior
changes the position of the microbe, so later behaviors see fresh values.
"""
invalidate_field_cache!(model::ABM) = (field_cache(model).id = 0; nothing)

@inline function _cached(compute::F, name::Symbol, microbe, model) where {F}
    cache = field_cache(model)
    cache.id == microbe.id || return compute()
    v = getfield(cache, name)
    v === nothing || return v
    v = compute()
    setfield!(cache, name, v)
    return v
end

function concentration(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:concentration, microbe, model) do
        concentration(chemoattractant(model))(microbe, model)::Float64
    end
end
function gradient(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:gradient, microbe, model) do
        gradient(chemoattractant(model))(microbe, model)::SVector{D,Float64}
    end
end
function time_derivative(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:time_derivative, microbe, model) do
        time_derivative(chemoattractant(model))(microbe, model)::Float64
    end
end
function chemoattractant_diffusivity(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    _cached(:chemoattractant_diffusivity, microbe, model) do
        chemoattractant_diffusivity(chemoattractant(model))(microbe, model)::Float64
    end
end

"""
    chemoattractant(model)
Returns the chemoattractant object from `model`.
"""
chemoattractant(model::ABM) = model.chemoattractant
"""
    concentration(model)
Returns the function `f` that defines the concentration field.
The returned function has signature `f(pos, model)` and returns a scalar.
"""
concentration(model::ABM) = concentration(chemoattractant(model))
"""
    gradient(model)
Returns the function `f` that defines the gradient of the concentration field.
The returned function has signature `f(pos, model)` and returns a `SVector`
with the same dimensionality as the microbe position `pos`.
"""
gradient(model::ABM) = gradient(chemoattractant(model))
"""
    time_derivative(model)
Returns the function `f` that defines the time derivative of the concentration field.
The returned function has signature `f(pos, model)` and returns a scalar.
"""
time_derivative(model::ABM) = time_derivative(chemoattractant(model))
"""
    chemoattractant_diffusivity(model)
Returns the thermal diffusivity of the chemoattractant compound.
"""
chemoattractant_diffusivity(model::ABM) = chemoattractant_diffusivity(chemoattractant(model))
concentration(c::AbstractChemoattractant) = c.concentration_field
gradient(c::AbstractChemoattractant) = c.concentration_gradient
time_derivative(c::AbstractChemoattractant) = c.concentration_ramp
chemoattractant_diffusivity(c::AbstractChemoattractant) = c.diffusivity

"""
    GenericChemoattractant{D} <: AbstractChemoattractant{D}
Type for a generic chemoattractant field.
Field, gradient and ramp default to 0 everywhere in the domain.
Diffusivity defaults to 608 μm²/s everywhere in the domain. (Do not use 0 here, it may mess up some calculations)
"""
struct GenericChemoattractant{D} <: AbstractChemoattractant{D}
    concentration_field::Function
    concentration_gradient::Function
    concentration_ramp::Function
    diffusivity::Function
end
function GenericChemoattractant{D}(;
    concentration_field = (::AbstractMicrobe, ::ABM) -> zero(Float64), # μM
    concentration_gradient = (::AbstractMicrobe, ::ABM) -> zero(SVector{D,Float64}), # μM/μm
    concentration_ramp = (::AbstractMicrobe, ::ABM) -> zero(Float64), # μM/s
    diffusivity = (::AbstractMicrobe, ::ABM) -> Float64(608), # μm²/s
) where D
    GenericChemoattractant{D}(
        concentration_field, concentration_gradient, concentration_ramp, diffusivity
    )
end
