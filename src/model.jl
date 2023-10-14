# extended from Agents.jl
export StandardABM, UnremovableABM, add_agent!, run!
# exported from Agents.jl without extensions
export add_agent_pos!, ContinuousSpace

"""
    UnremovableABM(MicrobeType, space, timestep; kwargs...)
Extension of the `Agents.UnremovableABM` method for microbe types.
Implementation of `AgentBasedModel` where agents can only be added but not removed.
See `Agents.AgentBasedModel` for detailed information on the keyword arguments.

**Arguments**
- `MicrobeType`: subtype of `AbstractMicrobe{D}`, with explicitly specified dimensionality `D`. A list of available options can be obtained by running `subtypes(AbstractMicrobe)`.
- `space`: a `ContinuousSpace{D}` with _the same_ dimensionality `D` as MicrobeType which specifies the spatial properties of the simulation domain.
- `timestep`: the integration timestep of the simulation.

**Keywords**
- `properties`: additional container of data to specify model-level properties. MicrobeAgents.jl includes a set of default properties (detailed at the end).
- `scheduler = Schedulers.fastest`
- `rng = Random.default_rng()`
- `spacing = minimum(extent)/20`
- `warn = true`

**Default `properties`**

When a model is created, a default set of properties is included in the model
(`MicrobeAgents.default_ABM_properties`):
```
DEFAULT_ABM_PROPERTIES = Dict(
    :t => 0, # counter for timekeeping
    :concentration_field => (pos,model) -> 0.0,
    :concentration_gradient => (pos,model) -> zero.(pos),
    :concentration_time_derivative => (pos,model) -> 0.0,
    # required by models of chemotaxis, default value is glutamate diffusivity
    :compound_diffusivity => 608.0,
    # model stepper, by default only keeps time
    :update! => tick!
)
```
By including these default properties, we make sure that all the chemotaxis models
will work even without extra user intervention.
All these properties can be overwritten by simply passing an equivalent key
to the `properties` dictionary when creating the model.
"""
function Agents.UnremovableABM(
    T::Type{A}, space::ContinuousSpace{D}, timestep::Real;
    scheduler = Schedulers.fastest,
    properties = Dict(),
    rng = Random.default_rng(),
    warn = true
) where {D,A<:AbstractMicrobe{D}}
    properties = Dict(
        DEFAULT_ABM_PROPERTIES...,
        properties...,
        :timestep => timestep
    )
    UnremovableABM(T, space; scheduler, properties, rng, warn)
end


"""
    StandardABM(MicrobeType, space, timestep; kwargs...)
Extension of the `Agents.StandardABM` method for microbe types.
Implementation of `AgentBasedModel` where agents can be added and removed at any time.
If agents removal is not required, it is recommended to use `UnremovableABM` for better performance.
See `Agents.AgentBasedModel` for detailed information on the keyword arguments.

**Arguments**
- `MicrobeType`: subtype of `AbstractMicrobe{D}`, with explicitly specified dimensionality `D`. A list of available options can be obtained by running `subtypes(AbstractMicrobe)`.
- `space`: a `ContinuousSpace{D}` with _the same_ dimensionality `D` as MicrobeType which specifies the spatial properties of the simulation domain.
- `timestep`: the integration timestep of the simulation.

**Keywords**
- `properties`: additional container of data to specify model-level properties. MicrobeAgents.jl includes a set of default properties (detailed at the end).
- `scheduler = Schedulers.fastest`
- `rng = Random.default_rng()`
- `warn = true`

**Default `properties`**

When a model is created, a default set of properties is included in the model
(`MicrobeAgents.default_ABM_properties`):
```
DEFAULT_ABM_PROPERTIES = Dict(
    :t => 0, # counter for timekeeping
    :concentration_field => (pos,model) -> 0.0,
    :concentration_gradient => (pos,model) -> zero.(pos),
    :concentration_time_derivative => (pos,model) -> 0.0,
    # required by models of chemotaxis, default value is glutamate diffusivity
    :compound_diffusivity => 608.0,
    # model stepper, by default only keeps time
    :update! => tick!
)
```
By including these default properties, we make sure that all the chemotaxis models
will work even without extra user intervention.
All these properties can be overwritten by simply passing an equivalent key
to the `properties` dictionary when creating the model.
"""
function Agents.StandardABM(
    T::Type{A}, space::ContinuousSpace{D}, timestep::Real;
    scheduler = Schedulers.fastest,
    properties = Dict(),
    rng = Random.default_rng(),
    warn = true,
) where {D,A<:AbstractMicrobe{D}}
    properties = Dict(
        DEFAULT_ABM_PROPERTIES...,
        properties...,
        :timestep => timestep
    )
    StandardABM(T, space; scheduler, properties, rng, warn)
end


# Microbe constructor generates a random velocity in non-reproducible way.
# When microbes are created internally these velocity must be generated
# reproducibly using the model rng, if a vel keyword is not specified.
# In this case, we want to initialize the microbe with a zero velocity
# and only then extract a random velocity with the correct rng.
# It is sufficient to extend this single method because it is
# the lowest level method to which all the others fall back.
"""
    add_agent!([pos,] [MicrobeType,] model; kwargs...)
MicrobeAgents extension of `Agents.add_agent!`.
Creates and adds a new microbe to `model`, using the constructor of the agent type
of the model.
If `model` accepts mixed agent types, then `MicrobeType` must be specified.
If not specified, `pos` will be assigned randomly in the model domain.

Keywords can be used to specify default values to pass to the microbe constructor,
otherwise default values from the constructor will be used.
If unspecified, a random velocity vector and a random speed are generated.
"""
function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{<:AbstractMicrobe{D}},
    model::AgentBasedModel,
    properties...;
    vel = nothing,
    speed = nothing,
    kwproperties...
) where {D}
    id = nextid(model)
    if !isempty(properties)
        microbe = A(id, pos, properties...)
    else
        microbe = A(; id, pos, vel = zero(SVector{D}), speed = 0.0, kwproperties...)
        microbe.vel = isnothing(vel) ? random_velocity(model) : vel
        microbe.speed = isnothing(speed) ? random_speed(microbe, model) : speed
    end
    add_agent_pos!(microbe, model)
end


"""
    run!(model, n=1; kwargs...)
MicrobeAgents extension of `Agents.run!`, which assumes the `agent_step!` to be
`microbe_step!`, and the `model_step!` to be `model.update!`.
Runs the model for a number of steps specified by `n`. If `n` is not specified, only 1 step is performed.
`n` can also be a function `n(model,s)::Bool` (where `s` is the current number of steps taken)
in which case the simulations stops when `n` returns true.

**Keywords**
- `when = true`: at which steps to perform collection/processing of agent data.
- `when_model = when`: at which steps to perform collection/processing of model data.
- `adata::Vector`: agent data to collect
- `mdata::Vector`: model data to collect
- `obtainer = identity`
- `agents_first = true`: whether agent stepping should be performed before model stepping
- `showprogress = false`
See `Agents.run!` and `Agents.step!` for detailed information.
"""
function Agents.run!(model::AgentBasedModel{S}, n = 1;
    when = true,
    when_model = when,
    adata = nothing,
    mdata = nothing,
    obtainer = identity,
    agents_first = true,
    showprogress = false
) where {S<:ContinuousSpace}
    run!(model, microbe_step!, model.update!, n;
        when, when_model,
        adata, mdata,
        obtainer, agents_first, showprogress
    )
end


"""
    tick!(model::AgentBasedModel)
Increase time count `model.t` by 1.
"""
tick!(model::AgentBasedModel) = (model.t += 1)
# extend function chaining
→(model::AgentBasedModel, f, g...) = (model.update! = →(model.update! → f, g...))
→(model::AgentBasedModel, f) = (model.update! = model.update! → f)

DEFAULT_ABM_PROPERTIES = Dict(
    :t => 0, # counter for timekeeping
    :concentration_field => (pos,model) -> 0.0,
    :concentration_gradient => (pos,model) -> zero.(pos),
    :concentration_time_derivative => (pos,model) -> 0.0,
    # required by models of chemotaxis, default value is glutamate diffusivity
    :compound_diffusivity => 608.0,
    # model stepper, by default only keeps time
    :update! => tick!
)
