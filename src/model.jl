"""
    StandardABM(MicrobeType, space, timestep; kwargs...)
Extension of the `Agents.StandardABM` method for microbe types.
Implementation of `AgentBasedModel` where agents can be added and removed at any time.
If agents removal is not required, it is recommended to use the
keyword argument `container = Vector` for better performance.
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
    :chemoattractant => GenericChemoattractant{D,Float64}()
    :affect! => chemotaxis!
)
```
By including these default properties, we make sure that all the chemotaxis models
will work even without extra user intervention.
All these properties can be overwritten by simply passing an equivalent key
to the `properties` dictionary when creating the model.
"""
function Agents.StandardABM(
    T::Type{A}, space::ContinuousSpace{D}, timestep::Real;
    agent_step! = microbe_step!,
    model_step! = _ -> nothing,
    container = Dict,
    scheduler = Schedulers.fastest,
    properties = Dict(),
    rng = Random.default_rng(),
    agents_first = true,
    warn = true,
) where {D,A<:AbstractMicrobe{D}}
    properties = (;
        make_default_abm_properties(D)...,
        properties...,
        timestep = timestep
    )
    StandardABM(T, space;
        agent_step!, model_step!, container,
        scheduler, properties, rng, agents_first, warn
    )
end


function Agents.add_agent!(
    pos::Agents.ValidPos,
    A::Type{<:AbstractMicrobe{D}},
    model::AgentBasedModel,
    properties...;
    vel = nothing,
    speed = nothing,
    kwproperties...
) where {D}
    @assert haskey(kwproperties, :motility) "Missing required keyword argument `motility`"
    N = get_motility_N(kwproperties[:motility])
    add_agent!(pos, A{N}, model, properties...; vel, speed, kwproperties...)
end

get_motility_N(m::Motility{N}) where {N} = N

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
    A::Type{<:AbstractMicrobe{D,N}},
    model::AgentBasedModel,
    properties...;
    vel = nothing,
    speed = nothing,
    kwproperties...
) where {D,N}
    @assert haskey(kwproperties, :motility) "Missing required keyword argument `motility`"
    id = Agents.nextid(model) # not public API!
    if !isempty(properties)
        microbe = A(id, pos, properties...)
    else
        microbe = A(; id, pos, vel = zero(SVector{D}), speed = 0.0, kwproperties...)
        microbe.vel = isnothing(vel) ? random_velocity(model) : vel
        microbe.speed = isnothing(speed) ? random_speed(microbe, model) : speed
    end
    Agents.add_agent_own_pos!(microbe, model) # not public API!
end

make_default_abm_properties(D) = Dict(
    :chemoattractant => GenericChemoattractant{D,Float64}(),
    :affect! => chemotaxis!
)
