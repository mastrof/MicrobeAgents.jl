export AbstractChemoattractant, GenericChemoattractant, chemoattractant
export concentration, gradient, time_derivative, chemoattractant_diffusivity

"""
    AbstractChemoattractant{D,T}
Abstract type for chemoattractants.
Requires dimensionality (`D`) and number type (`T`, typically `Float64`).

The interface is defined by four core functions that operate on `AgentBasedModel`s:
- `chemoattractant`: returns the chemoattractant object
- `concentration`: returns the function for the concentration field
- `gradient`: returns the function for the concentration gradient
- `time_derivative`: returns the function for the concentration ramp
- `chemoattractant_diffusivity`: returns the thermal diffusivity of the chemoattractant
"""
abstract type AbstractChemoattractant{D,T} end

@inline function concentration(pos::SVector{D,T}, model::ABM) where {D,T}
    concentration(chemoattractant(model))(pos, model)::T
end
@inline function gradient(pos::SVector{D,T}, model::ABM) where {D,T}
    gradient(chemoattractant(model))(pos, model)::SVector{D,T}
end
@inline function time_derivative(pos::SVector{D,T}, model::ABM) where {D,T}
    time_derivative(chemoattractant(model))(pos, model)::T
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
gradient(model::ABM) = concentration(chemoattractant(model))
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
    GenericChemoattractant{D,T} <: AbstractChemoattractant{D,T}
Type for a generic chemoattractant field.
Concentration field, gradient and time derivative default to zero values.
Diffusivity defaults to 608 μm^2/s.
"""
@kwdef struct GenericChemoattractant{D,T} <: AbstractChemoattractant{D,T}
    concentration_field::Function = (pos, model) -> zero(T)
    concentration_gradient::Function = (pos, model) -> zero(SVector{D,T})
    concentration_ramp::Function = (pos, model) -> zero(T)
    diffusivity::T = T(608.0) # μm²/s
end
