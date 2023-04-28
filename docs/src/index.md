# MicrobeAgents.jl

MicrobeAgents.jl is a Julia framework for agent-based simulations of bacterial
behavior, built on the amazing [Agents.jl](https://github.com/JuliaDynamics/Agents.jl).

## Features
- Runs in 1, 2 and 3 spatial dimensions.
- Provides various motility patterns (Run-Tumble, Run-Reverse, Run-Reverse-Flick), all with customizable speed and turn angle distributions.
- Various models of bacterial chemotaxis (Brown & Berg, PNAS 1974; Celani & Vergassola, PNAS 2010; Xie et al, Biophys J 2014; Brumley et al, PNAS 2019).
- Fast analysis routines for common quantities of interest (run statistics, MSD, autocorrelation functions, drift velocity).

## Limitations (some may be temporary, others may be not)
- Only continuous space models are supported
- Reorientations are always assumed to be instantaneous; this approximation is really only reasonable when the integration timestep is above 50ms.
- Integration timestep also sets the sensory integration timescale in chemotactic models.

## What this package is not good for
Although, in principle, you can add arbitrary layers of complexity on top the provided interface, there are a few things for which this package is not a recommended choice.
- Hydrodynamic interactions.
- Atomistic representation of biochemical pathways.

## Contribute
The package is still in an early stage of intense development.
If you would like to have support for your favorite model of chemotaxis, or need some specific features to be implemented, please open an issue. I'll try to satisfy as many requests as possible.

If you would like to take a more active part in the development, please consider contacting me directly at rfoffi@ethz.ch.

## Citation
If you use this package in work that leads to a publication, please cite the GitHub repository:
```
@misc{Foffi2023,
    author = {Foffi, R.},
    title = {MicrobeAgents.jl},
    year = {2023},
    publisher = {GitHub},
    journal = {GitHub repository},
    howpublished = {\url{https://github.com/mastrof/MicrobeAgents.jl}}
}
```

## Acknowledgements
This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 955910.