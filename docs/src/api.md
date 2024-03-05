# API

## [Microbes](@id Microbes)
```@docs
AbstractMicrobe
Microbe
BrownBerg
Brumley
Celani
Xie
```

```@docs
position
direction
speed
velocity
motilepattern
turnrate
rotational_diffusivity
radius
state
```

## [Motility](@id Motility)
```@docs
AbstractMotility
MotilityOneStep
MotilityTwoStep
RunTumble
RunReverse
RunReverseFlick
```

## [Utils](@id Utils)
```@docs
distance
distancevector
random_speed
random_velocity
Analysis.adf_to_matrix
Analysis.adf_to_vectors
Analysis.unfold
Analysis.unfold!
```

## [Data analysis](@id Analysis)
```@docs
Analysis.detect_turns!
Analysis.run_durations
Analysis.driftvelocity_point!
Analysis.driftvelocity_direction!
msd
acf
```
