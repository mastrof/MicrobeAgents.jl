export ABM, add_agent!, add_agent_pos!, run!

"""
    tick!(model::ABM)
Increase time count `model.t` by 1.
"""
tick!(model::ABM) = (model.t += 1)
# extend function chaining
→(model::ABM, f, g...) = (model.update! = →(model.update! → f, g...))
→(model::ABM, f) = (model.update! = model.update! → f)

default_ABM_properties = Dict(
    :t => 0, # counter for timekeeping
    :concentration_field => (pos,model) -> 0.0,
    :concentration_gradient => (pos,model) -> zero.(pos),
    :concentration_time_derivative => (pos,model) -> 0.0,
    # required by models of chemotaxis, default value is glutamate diffusivity
    :compound_diffusivity => 608.0,
    # model stepper, by default only keeps time
    :update! => tick!
)

function Agents.AgentBasedModel(
    T::Type{A},
    extent::NTuple{D,<:Real}, timestep::Real;
    periodic = true,
    scheduler::F = Schedulers.fastest,
    properties::P = Dict(),
    rng::R = Random.default_rng(),
    warn = true,
) where {D,A<:AbstractMicrobe{D},F,P,R<:AbstractRNG}
    space = ContinuousSpace(extent; periodic)
    properties = Dict(
        default_ABM_properties...,
        properties...,
        :timestep => timestep
    )
    ABM(T, space; scheduler, properties, rng, warn)
end

function Agents.AgentBasedModel(microbe::AbstractMicrobe, args...; kwargs...)
    return ABM(typeof(agent), args...; kwargs...)
end

# Microbe constructor generates a random velocity in non-reproducible way.
# When microbes are created internally these velocity must be generated
# reproducibly using the model rng, if a vel keyword is not specified.
# In this case, we want to initialize the microbe with a zero velocity
# and only then extract a random velocity with the correct rng.
# It is sufficient to extend this single method because it is
# the lowest level method to which all the others fall back.
function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{<:AbstractMicrobe{D}},
    model::ABM,
    properties...;
    vel = nothing,
    speed = nothing,
    kwargs...
) where D
    id = nextid(model)
    microbe = A(id, pos, properties...; vel=ntuple(zero,D), speed=0, kwargs...)
    microbe.vel = isnothing(vel) ? rand_vel(model.rng, D) : vel
    microbe.speed = isnothing(speed) ? rand_speed(model.rng, microbe.motility) : speed
    add_agent_pos!(microbe, model)
end

function Agents.run!(model::ABM{S,A,F,P,R}, n = 1;
    when = true,
    when_model = when,
    adata = nothing,
    mdata = nothing,
    obtainer = identity,
    agents_first = true,
    showprogress = false
) where {S<:ContinuousSpace,A<:AbstractMicrobe,F,P,R<:AbstractRNG}
    run!(model, microbe_step!, model.update!, n;
        when, when_model,
        adata, mdata,
        obtainer, agents_first, showprogress
    )
end