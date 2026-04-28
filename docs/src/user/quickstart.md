# Quickstart

Use this page if you want the fastest path from clone to a successful local
run.

This page is for first-time users who do not need GRAM or SPICE on day one.

This default quickstart baseline is the Earth vacuum path:
`NoAtmosphereModel()` plus `SimpleEphemeridesModel()`.

Shortest successful path:

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Assets & Modes](../assets.md)
- [First Simulation](first_simulation.md)
- [Recipes](recipes.md)

## What this run gives you

The no-GRAM example is the shortest repository-owned path to a working local
simulation. It does not require:

- a local GRAM installation
- SPICE kernels
- surrogate atmosphere grids

## If you prefer the CLI

The equivalent command-line path is:

```bash
./bin/spaceagora run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

## When to stop using this page

Move on as soon as one of these becomes true:

- you need higher-fidelity assets or licensed data
- you want to run a fuller scenario than the baseline no-GRAM example
- you want the verification study rather than a first smoke run
