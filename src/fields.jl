export AbstractChemoattractant, chemoattractant,
    concentration, gradient, time_derivative, chemoattractant_diffusivity,
    GenericChemoattractant

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

chemoattractant(model::ABM) = model.chemoattractant
concentration(model::ABM) = concentration(chemoattractant(model))
gradient(model::ABM) = concentration(chemoattractant(model))
time_derivative(model::ABM) = time_derivative(chemoattractant(model))
chemoattractant_diffusivity(model::ABM) = chemoattractant_diffusivity(chemoattractant(model))
concentration(c::AbstractChemoattractant) = c.concentration_field
gradient(c::AbstractChemoattractant) = c.concentration_gradient
time_derivative(c::AbstractChemoattractant) = c.concentration_ramp
chemoattractant_diffusivity(c::AbstractChemoattractant) = c.diffusivity

@kwdef struct GenericChemoattractant{D,T} <: AbstractChemoattractant{D,T}
    concentration_field::Function = (pos, model) -> zero(T)
    concentration_gradient::Function = (pos, model) -> zero(SVector{D,T})
    concentration_ramp::Function = (pos, model) -> zero(T)
    diffusivity::T = T(608.0) # μm²/s
end
