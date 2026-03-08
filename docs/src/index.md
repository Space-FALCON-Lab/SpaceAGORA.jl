# SpaceAGORA Documentation

SpaceAGORA is a Julia toolkit for spacecraft orbit and atmospheric-flight
simulation. This site is organized for three audiences:

- users who want a first successful run
- users who need command, API, or extension reference material
- maintainers working on package surface, documentation, and repository rules

## Start in under a minute

Use this path if you want the fastest successful run on a machine that does not
have GRAM or SPICE configured yet.

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=.AGORA -e 'using Pkg; Pkg.instantiate()'
julia --project=.AGORA examples/AGORA_Earth_NoGRAM.jl
```

## Choose your path

### No-GRAM first run

Use this route if you want the shortest path to a working local simulation.

- [Quickstart](user/quickstart.md)
- [Assets & Modes](assets.md)
- [First Simulation](user/first_simulation.md)

### High-fidelity GRAM/SPICE

Use this route if you need higher-fidelity atmosphere, ephemerides, or frame
data.

- [Installation & Environment](user/installation_environment.md)
- [Assets & Modes](assets.md)
- [Verification Study](user/verification_study.md)

## Documentation tracks

### User Guide

Start here if your goal is to run something, understand the main workflows, or
find a concrete next step.

- [Start Here](getting_started.md)
- [Quickstart](user/quickstart.md)
- [Recipes](user/recipes.md)

### Reference

Use the reference pages when you already know what you want to do and need the
exact command surface, extension hooks, or generated API details.

- [CLI](cli.md)
- [Public API](generated/public_api.md)
- [Extensibility](extensibility.md)

### Maintainer Guide

Use the maintainer section when you are changing docs, public API shape, repo
structure, or CI-backed documentation rules.

- [Maintainer Overview](maintainer/index.md)
- [Documentation Policy](documentation_policy.md)
- [Contracts](contracts.md)
