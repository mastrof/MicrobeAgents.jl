export ABM, add_agent!, add_agent_pos!, run!

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
"""
    tick!(model::ABM)
Increase time count `model.t` by 1.
"""
tick!(model::ABM) = model.t += 1
# extend function chaining
→(model::ABM, f, g...) = (model.update! = →(model.update! → f, g...))
→(model::ABM, f) = (model.update! = model.update! → f)

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