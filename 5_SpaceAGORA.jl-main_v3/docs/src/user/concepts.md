# Concepts

Use this page when you want the mental model behind the main workflows without
digging into maintainer-only architecture documents.

This page is for users who can already run commands and now want to understand
which route, mode, or surface is appropriate for their work.

Shortest successful command:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

What to read next:

- [First Simulation](first_simulation.md)
- [Assets & Modes](../assets.md)
- [Public API](../generated/public_api.md)

## Two operating modes

### No-GRAM mode

Use no-GRAM mode when you want:

- onboarding
- smoke runs
- lightweight studies
- the shortest route to a successful run on a new machine

### Higher-fidelity mode

Use the higher-fidelity path when you need:

- mission-quality atmosphere and winds
- SPICE-backed geometry and frame data
- other machine-local or licensed assets

## Three common entry surfaces

### Example scripts

Example scripts are the quickest way to run a repository-owned scenario.

### CLI

The CLI is the stable packaged command surface for examples, verification,
benchmarks, and asset inspection.

### Root package API

For scripts and extension work, depend on the exported `SpaceAGORA` surface.
Treat the generated Public API page as the supported list of exports.

## Outputs stay local

SpaceAGORA regenerates reports, docs builds, and similar artifacts into ignored
paths such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.

## Performance knob: `isolate_state`

`SpaceAGORA.run_simulation(...; isolate_state=true)` deep-copies the
configuration by default so repeated or concurrent runs do not alias mutable
state.

Treat `isolate_state=false` as an expert-only throughput lever after you have
measured the setup cost and confirmed that the configuration instance is not
shared across overlapping or state-mutating runs.
