## Creating a microbe
Microbes are represented by subtypes of the `AbstractMicrobe` type, which is in turn a subtype of `AbstractAgent` introduced by Agents.jl
```@docs
AbstractMicrobe
```

MicrobeAgents provides different `AbstractMicrobe` subtypes representing different models of bacterial behavior from the literature.

A basic type, which is typically sufficient for simple motility simulations and does not include chemotaxis, is the `Microbe` type.
```@docs
Microbe
```

The dimensionality of `Microbe` *must* always be specified on creation. All the fields are instead optional, and if not specified will be assigned default values.

To create a `Microbe` living in a 1-dimensional space, with run-tumble motility and average turn rate ``\nu=1\;s^{-1}``, it is therefore sufficient to run
```
Microbe{1}()
```
Similarly, for two and three dimensions:
```
Microbe{2}()
Microbe{3}()
```

Any custom parameter can be set via kwargs:
```
Microbe{3}(
    turn_rate = 0.6,
    rotational_diffusivity = 0.1
)
```


All the other subtypes of `AbstractMicrobe` work in a similar way, although they will have distinct default values and extra fields.

```@docs
BrownBerg
Brumley
Celani
Xie
```


## Creating a model
MicrobeAgents.jl exploits the `AgentBasedModel` interface from Agents.jl.
While the standard Agents.jl syntax will always work, it is typically more convenient to use the method extensions provided by MicrobeAgents.jl, which also includes some default parameters required by the simulations.
Both `StandardABM` and `UnremovableABM` are supported.


Whenever removal of microbes during the simulation is not needed, `UnremovableABM` is the recommended choice.