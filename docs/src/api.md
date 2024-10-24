# API

## [Microbes](@id Microbes)
```@docs
position
direction
speed
velocity
motilepattern
rotational_diffusivity
radius
state
motilestate
states
transition_weights
duration
angle
distance
distancevector
random_speed
random_velocity
```

## [Chemoattractants](@id Chemoattractants)
```@docs
AbstractChemoattractant
GenericChemoattractant
chemoattractant
concentration
gradient
time_derivative
chemoattractant_diffusivity
```

## [Utils](@id Utils)
```@docs
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
