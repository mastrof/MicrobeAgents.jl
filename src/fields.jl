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

function concentration(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    concentration(chemoattractant(model))(microbe, model)::Float64
end
function gradient(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    gradient(chemoattractant(model))(microbe, model)::SVector{D,Float64}
end
function time_derivative(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    time_derivative(chemoattractant(model))(microbe, model)::Float64
end
function chemoattractant_diffusivity(microbe::AbstractMicrobe{D,N}, model::ABM) where {D,N}
    chemoattractant_diffusivity(chemoattractant(model))(microbe, model)::Float64
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
    GenericChemoattractant{D} <: AbstractChemoattractant{D}
Type for a generic chemoattractant field.
Field, gradient and ramp default to 0 everywhere in the domain.
Diffusivity defaults to 608 μm²/s everywhere in the domain. (Do not use 0 here, it may mess up some calculations)
"""
@kwdef struct GenericChemoattractant{D} <: AbstractChemoattractant{D}
    concentration_field::Function = (::AbstractMicrobe, ::ABM) -> zero(Float64) # μM
    concentration_gradient::Function = (::AbstractMicrobe, ::ABM) -> zero(SVector{D,Float64}) # μM/μm
    concentration_ramp::Function = (::AbstractMicrobe, ::ABM) -> zero(Float64) # μM/s
    diffusivity::Function = (::AbstractMicrobe, ::ABM) -> Float64(608) # μm²/s
end
