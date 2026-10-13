# MicrobeAgents.jl

MicrobeAgents.jl is a Julia framework for agent-based simulations of bacterial
motility and chemotaxis, built on the amazing [Agents.jl](https://github.com/JuliaDynamics/Agents.jl).

## Features
- Runs in 1, 2 and 3 spatial dimensions.
- Provides base motility patterns (Run-Tumble, Run-Reverse, Run-Reverse-Flick, Run-Stop), all with customizable speed and turn angle distributions, and allows user definition of new arbitrary patterns.
- Includes various models of bacterial chemotaxis (Brown & Berg, PNAS 1974; Celani & Vergassola, PNAS 2010; Xie et al, Biophys J 2014; Brumley et al, PNAS 2019; Son, Menolascina & Stocker, PNAS 2016).
- Modular and customizable interface for microbial behaviors
- Different options to represent physically structured environments
- Fast analysis routines for common quantities of interest (MSD, autocorrelation functions, drift velocity).

## Limitations (some may be temporary, others may be not)
- Only continuous space models are supported
- Integration timestep also sets the sensory integration timescale in chemotactic models.
- The library adopts fundamentally a stochastic kinetic approach. Full integration of dynamic (i.e., force- and torque-based) equations of motion is not currently in scope.

## Contribute
If you would like to have support for your favorite model of chemotaxis, or need some specific features to be implemented, please open an issue.

If you would like to take a more active part in the development, please consider contacting me directly at rfoffi@ethz.ch.

## Citation
If you use this package in work that leads to a publication, please cite the
[Zenodo record](https://doi.org/10.5281/zenodo.14786182) of the software
(the DOI always resolves to the latest release):

```@eval
using MicrobeAgents, Dates, Markdown
v = pkgversion(MicrobeAgents)
fence = "`"^3
Markdown.parse("""
$(fence)
@software{Foffi_MicrobeAgents,
    author = {Foffi, Riccardo},
    title = {MicrobeAgents.jl},
    version = {v$(v)},
    year = {$(Dates.year(Dates.today()))},
    publisher = {Zenodo},
    doi = {10.5281/zenodo.14786182},
    url = {https://doi.org/10.5281/zenodo.14786182}
}
$(fence)
""")
```

## Acknowledgements
This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 955910.
