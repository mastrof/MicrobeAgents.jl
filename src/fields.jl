export AbstractConcentrationField, field, concentration, gradient, time_derivative,
    GenericConcentrationField,
    SteadyDiffusiveField, background, intensity, source_size, origin,
    SteadyExpField, decay_length

abstract type AbstractConcentrationField{D} end

field(model::ABM) = model.field
concentration(model::ABM) = concentration(field(model))
gradient(model::ABM) = gradient(field(model))
time_derivative(model::ABM) = time_derivative(field(model))
concentration(f::AbstractConcentrationField) = f.concentration_field
gradient(f::AbstractConcentrationField) = f.concentration_gradient
time_derivative(f::AbstractConcentrationField) = f.concentration_ramp

@kwdef struct GenericConcentrationField{D} <: AbstractConcentrationField{D}
    concentration_field = (pos, model) -> 0.0
    concentration_gradient = (pos, model) -> zero(SVector{D,Float64})
    concentration_ramp = (pos, model) -> 0.0
end


@kwdef struct SteadyDiffusiveField{D} <: AbstractConcentrationField{D}
    C0::Float64
    C1::Float64
    R::Float64
    origin::SVector{D,Float64} = zero(SVector{D,Float64})
    concentration_field = field_diffusive
    concentration_gradient = gradient_diffusive
    concentration_ramp = (pos, model) -> 0.0
end
background(f::SteadyDiffusiveField) = f.C0
intensity(f::SteadyDiffusiveField) = f.C1
source_size(f::SteadyDiffusiveField) = f.R
origin(f::SteadyDiffusiveField) = f.origin

function field_diffusive(pos::SVector{D}, model) where D
    F = field(model)
    C0 = background(F)
    C1 = intensity(F)
    R = source_size(F)
    P = origin(F)
    r = euclidean_distance(pos, P, model)
    return C0 + C1*R/r
end

function gradient_diffusive(pos::SVector{D}, model) where D
    F = field(model)
    C1 = intensity(F)
    R = source_size(F)
    P = origin(F)
    rvec = distancevector(P, pos, model)
    r2 = dot(rvec, rvec)
    r3 = r2 * sqrt(r2)
    return SVector{D}(-C1*R*r/r3 for r in rvec)
end


@kwdef struct SteadyExpField{D} <: AbstractConcentrationField{D}
    C0::Float64
    C1::Float64
    λ::Float64
    R::Float64
    origin::SVector{D,Float64} = zero(SVector{D,Float64})
    concentration_field = field_exp
    concentration_gradient = gradient_exp
    concentration_ramp = (pos, model) -> 0.0
end
background(f::SteadyExpField) = f.C0
intensity(f::SteadyExpField) = f.C1
source_size(f::SteadyExpField) = f.R
origin(f::SteadyExpField) = f.origin
decay_length(f::SteadyExpField) = f.λ

function field_exp(pos::SVector{D}, model) where D
    F = field(model)
    C0 = background(F)
    C1 = intensity(F)
    R = source_size(F)
    P = origin(F)
    λ = decay_length(F)
    r = euclidean_distance(pos, P, model)
    return C0 + C1*R/r * exp(-(r-R)/λ)
end

function gradient_exp(pos::SVector{D}, model) where D
    F = field(model)
    C1 = intensity(F)
    R = source_size(F)
    P = origin(F)
    λ = decay_length(F)
    rvec = distancevector(P, pos, model)
    r2 = dot(rvec, rvec)
    r = sqrt(r2)
    r3 = r*r2
    SVector{D}(-(λ+r)*r*exp(-(r-R)/λ)/(λ*r3) for r in rvec)
end
