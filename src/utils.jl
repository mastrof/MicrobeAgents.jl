export random_velocity, random_speed

"""
    random_velocity(model)
Generate a random velocity vector with unit norm respecting the
dimensionality of `model`.
"""
function random_velocity(model::AgentBasedModel{<:ContinuousSpace{D}}) where D
    random_velocity(abmrng(model), D)
end
function random_velocity(rng::AbstractRNG, D::Int)
    v = randn(rng, SVector{D})
    mag = sqrt(dot(v,v))
    return v ./ mag
end

"""
    random_speed(microbe, model)
Generate a random speed from the motile pattern of `microbe`.
"""
function random_speed(microbe::AbstractMicrobe, model::AgentBasedModel)
    rand(abmrng(model), speed(motilepattern(microbe)))
end
