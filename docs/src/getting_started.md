# Start Here

Use this page if you are new to SpaceAGORA and want the shortest route to the
right next document.

This page is for first-time users choosing between a no-GRAM local run, a
higher-fidelity setup, or a packaged study workflow.

The root environment is the canonical committed execution environment for
the repository-owned quickstart path.

The default minimal onboarding path is the no-GRAM Earth vacuum baseline:
`NoAtmosphereModel()` with `SimpleEphemeridesModel()`.

Minimal API check:

```jldoctest
julia> using SpaceAGORA

julia> env = make_no_gram_environment();

julia> env.density_model isa NoAtmosphereModel
true
```

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Quickstart](user/quickstart.md)
- [Assets & Modes](assets.md)
- [First Simulation](user/first_simulation.md)

## Pick the next page

If you want the fastest local run:

- go to [Quickstart](user/quickstart.md)

If you need setup details for the repository environments:

- go to [Installation & Environment](user/installation_environment.md)

If you need to know whether your machine has the right data or licensed assets:

- go to [Assets & Modes](assets.md)

If you want to run a fuller scenario after the first smoke path:

- go to [First Simulation](user/first_simulation.md)

If you want the telemetry verification workflow:

- go to [Verification Study](user/verification_study.md)

If you want copy-paste commands for common tasks:

- go to [Recipes](user/recipes.md)

## Before you dive deeper

SpaceAGORA supports two common working modes:

- a no-GRAM path for onboarding, smoke runs, and lightweight studies
- a higher-fidelity path that expects GRAM/SPICE assets on the local machine

You do not need the higher-fidelity path to get a first successful run.

## If you are maintaining the docs or repo structure

The user guide is not the right place for documentation policy or architecture
rules. Use the [Maintainer Overview](maintainer/index.md) for that material.
