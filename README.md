# MicrobeAgents.jl

[![DOI](https://zenodo.org/badge/587501553.svg)](https://doi.org/10.5281/zenodo.14786182)
[![Build Status](https://github.com/mastrof/MicrobeAgents.jl/workflows/CI/badge.svg)](https://github.com/mastrof/MicrobeAgents.jl/actions)
[![codecov](https://codecov.io/gh/mastrof/MicrobeAgents.jl/branch/main/graphs/badge.svg)](https://codecov.io/gh/mastrof/MicrobeAgents.jl)
[![Documentation, stable](https://img.shields.io/badge/docs-latest-blue.svg)](https://mastrof.github.io/MicrobeAgents.jl/dev/)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

MicrobeAgents.jl is a Julia framework for agent-based
simulations of microbial behavior (especially bacteria), built on
the amazing [Agents.jl](https://github.com/JuliaDynamics/Agents.jl).

Contributions, requests and suggestions are more than welcome.

## Main features
- Multiple swimming strategies (run-tumble, run-reverse, run-reverse-flick, run-stop) with tunable parameters, and possibility to define custom strategies
- Classical and modern models of chemotaxis (BrownBerg, Celani, Xie, Brumley, SonMenolascina)
- Support for arbitrary concentration fields, both numerical and analytical
- Modular behavioral traits (chemotaxis, chemokinesis...), also freely composable with custom user-defined ones
- Compatible with DifferentialEquations.jl for parallel integration of external fields and bacterial behavior
- Motility in complex environments through Agents.Pathfinding
- Analysis routines for standard quantities of interest (MSD, autocorrelation functions)

## Contribute
If you want to point out a bug, request some features or simply ask for info,
please don't hesitate to open an issue!

If you are interested in taking on a more active part in the development,
consider contacting me directly at rfoffi@ethz.ch.
I'll be happy to have a chat!
