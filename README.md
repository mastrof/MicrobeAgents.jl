# MicrobeAgents.jl

[![Build Status](https://github.com/mastrof/MicrobeAgents.jl/workflows/CI/badge.svg)](https://github.com/mastrof/MicrobeAgents.jl/actions)
[![codecov](https://codecov.io/gh/mastrof/MicrobeAgents.jl/branch/main/graphs/badge.svg)](https://codecov.io/gh/mastrof/MicrobeAgents.jl)
[![Documentation, stable](https://img.shields.io/badge/docs-latest-blue.svg)](https://mastrof.github.io/MicrobeAgents.jl/dev/)

MicrobeAgents.jl (previously Bactos.jl) is a Julia framework for agent-based
simulations of microbial behavior (especially bacteria), built on
the amazing [Agents.jl](https://github.com/JuliaDynamics/Agents.jl).
MicrobeAgents.jl extends and re-exports a minimal set of Agents.jl
functions and structures to be used as a stand-alone package, but it is
recommended to use it alongside Agents.jl for extra niceties.

The package is still at an early stage of intense development.
Contributions, requests and suggestions are more than welcome.

**To use the latest development version of MicrobeAgents it's required to manually install
the latest version of Agents first (`]add Agents#main`).**

**Documentation deployment is currently broken due to dependency on unreleased Agents.jl version; you can build it locally with the following code (requires the LiveServer.jl package):**
```
git clone https://github.com/mastrof/MicrobeAgents.jl.git
cd MicrobeAgents.jl
julia --project=docs/ -e 'import Pkg; Pkg.instantiate(); Pkg.add(name="Agents", rev="main"); include("docs/make.jl")'
julia -e 'using LiveServer; serve(dir="docs/build")'
```

## Main features
- Multiple swimming strategies (run-tumble, run-reverse, run-reverse-flick) with tunable parameters, and possibility to implement your own with minimal effort
- Classical and modern models of chemotaxis (Berg-Purcell, Celani-Vergassola, Xie, Brumley)
- Support for arbitrary concentration fields, both numerical and analytical
- Compatible with DifferentialEquations.jl for parallel integration of external fields and bacterial behavior
- Motility in complex environments through Agents.Pathfinding
- Analysis routines for standard quantities of interest (MSD, autocorrelation functions, drift velocity)

## Contribute
If you want to point out a bug, request some features or simply ask for info,
please don't hesitate to open an issue!

If you are interested in taking on a more active part in the development,
consider contacting me directly at rfoffi@ethz.ch.
I'll be happy to have a chat!
