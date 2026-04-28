# SpaceAGORA.jl

SpaceAGORA.jl is a Julia toolkit for spacecraft orbit and atmospheric-flight
simulation, including aerobraking, entry, and aerocapture workflows. The main
user-facing surfaces are:

- the root `SpaceAGORA` Julia API
- the repository CLI wrapper at `./bin/spaceagora`
- runnable example scripts under `examples/`

## Installation

Use the repository root environment for examples, tests, and normal local runs:

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
git submodule update --init --recursive --remote
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

That gives you the baseline open-data onboarding path immediately. You do not
need GRAM or SPICE to run the first quickstart example.

Helpful first checks:

```bash
./bin/spaceagora assets check
./bin/spaceagora assets setup-open
```

`assets check` reports which optional local assets are present.
`assets setup-open` summarizes the baseline repo-only path and which licensed
external assets remain user-provided.

## Asset model

SpaceAGORA supports two common setup modes:

### 1. Baseline no-GRAM mode

This mode is fully repository-local and is the recommended first run.

- uses built-in planets such as `Earth()`, `Mars()`, and `Venus()`
- uses fallback atmospheres such as `NoAtmosphereModel()` or
  `ExponentialAtmosphereModel(planet)`
- uses `SimpleEphemeridesModel()` instead of SPICE-backed ephemerides

### 2. GRAM/SPICE-backed mode

This mode is for higher-fidelity atmosphere and ephemeris workflows.

- requires a usable GRAM installation under
  `data/GRAMSuite.jl/GRAM Suite 2.0`
- requires local SPICE kernels under
  `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE`
- loads the `GRAMSuite` package extension automatically when `GRAMSuite` is
  importable

Need the exact process for pulling the `GRAMSuite.jl` submodule and copying the
official NASA GRAM Suite folders into the expected repo path after access is
approved? See
[docs/src/user/gramsuite_setup.md](docs/src/user/gramsuite_setup.md). It is a
step-by-step local setup guide for licensed GRAM assets.

If the native GRAM shared library for your operating system is missing, build
it with:

```bash
./scripts/ensure_gram_native.sh
```

If the build metadata came from a different machine or checkout path, force a
clean rebuild with:

```bash
./scripts/ensure_gram_native.sh --clean
```

## Quick Start

The fastest first run is:

```bash
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

Or through the CLI:

```bash
./bin/spaceagora run --example=AGORA_Basic_Quickstart.jl
```

That example is intentionally simple:

- Earth orbit
- no atmosphere drag
- no SPICE dependency
- no GRAM dependency
- results written under `output/`

If you only want to verify the end-to-end path with a short run, use:

```bash
./bin/spaceagora run --example=AGORA_Basic_Quickstart.jl --smoke
```

## Basic Example

The quickstart script lives at
`examples/AGORA_Basic_Quickstart.jl` and is the recommended first code sample
to read. It demonstrates the minimum pieces needed to run a simulation:

1. activate the shared example helpers
2. choose a no-GRAM planet model
3. build a simple spacecraft
4. assemble a `SimulationConfiguration`
5. run the simulation and save a CSV report

Key ideas in that example:

- `make_no_gram_planet(:earth)` avoids SPICE and uses built-in Earth constants
- `NoAtmosphereModel()` gives a vacuum baseline so installation is easier to
  validate
- `SimpleEphemeridesModel()` keeps the frame/ephemeris path open-data and
  lightweight
- `make_three_body_spacecraft(...)` builds a simple bus-plus-panels geometry
- `make_example_config(...)` packages common example defaults into a single
  configuration object

## GRAM Access

Once GRAM assets are available locally, the smallest GRAM-backed example is:

```bash
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

That example adds the minimum GRAM-specific pieces on top of the baseline path:

1. `setup_gram_example!()` loads the vendored `GRAMSuite` project if needed
2. `Earth("", SPICE_PATH)` uses the SPICE-backed planet constructor
3. `GRAMAtmosphereModel(planet_name="earth")` switches the atmosphere source
   from the built-in fallback model to GRAM

Before running it, verify:

- `data/GRAMSuite.jl/GRAM Suite 2.0` exists
- `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE` exists
- the platform-native GRAM shared library exists, or can be built with
  `./scripts/ensure_gram_native.sh`

Useful commands:

```bash
./bin/spaceagora assets check
./scripts/ensure_gram_native.sh
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

If you want fuller mission examples after that, start with:

- `examples/AGORA_Vex.jl`
- `examples/AGORA_Odyssey.jl`

Both of those use the same GRAM access pattern as the basic GRAM example, but
add mission-specific guidance, maneuvers, and higher-fidelity dynamics.

## CLI and docs

- Docs site: [space-falcon-lab.github.io/SpaceAGORA.jl](https://space-falcon-lab.github.io/SpaceAGORA.jl)
- CLI wrapper: `./bin/spaceagora`
- Run an example: `./bin/spaceagora run --example=AGORA_Basic_Quickstart.jl`
- Asset check: `./bin/spaceagora assets check`
- Verification study: `./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1`

## For contributors

- User docs live in `docs/src/` and are built with `julia --project=docs docs/make.jl`
- Generated API pages come from `docs/public_api_symbols.jl`
- Architecture and quality references live under the docs Maintainer Guide

Generated reports, docs builds, and other local outputs stay in ignored paths
such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.
