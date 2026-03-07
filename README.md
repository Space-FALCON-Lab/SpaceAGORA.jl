# SpaceAGORA.jl

## What SpaceAGORA is

SpaceAGORA.jl is a Julia toolkit for high-fidelity spacecraft orbit and atmospheric-flight simulation, including aerobraking, entry, and aerocapture scenarios. It supports typed simulation configuration, SciML-based execution through `run_simulation`, optional GRAM-backed atmospheric modeling, SPICE-backed ephemerides and frames, telemetry verification studies, and a package-owned CLI.

Generated API documentation and architecture references are published at [space-falcon-lab.github.io/SpaceAGORA.jl](https://space-falcon-lab.github.io/SpaceAGORA.jl).

The stable package contract is the root `SpaceAGORA` API, matching the repository's public API policy. For most users, the supported root entrypoints are `run_simulation`, `simulation_engine_config_from_env`, `run_cli`, `make_no_gram_planet`, `make_no_gram_density_model`, and `make_no_gram_environment`.

Internal modules such as `SimulationModel`, `SimulationEngine`, `ParallelProfiles`, `TelemetryVerification`, `RuntimeServices`, and `SpaceAGORACLI` are implementation detail namespaces, not the stable user-facing contract. If a symbol is not re-exported from `SpaceAGORA` and documented on the generated Public API page, treat it as internal.

## Quickstart with `.AGORA`

`/.AGORA` is the canonical committed execution environment for examples, benchmarks, tests, and CI. It is tracked in the repository, so there is no bootstrap step that copies the root `Project.toml` or `Manifest.toml` at runtime.

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=.AGORA -e 'using Pkg; Pkg.instantiate()'
```

Operational commands in this repository should continue to use `--project=.AGORA` unless that contract is deliberately changed across the repo.

## No-GRAM first run

SpaceAGORA supports a first-class no-GRAM onboarding path. This mode uses:

- `NoAtmosphereModel()` or `ExponentialAtmosphereModel(planet)`
- `SimpleEphemeridesModel()`
- planet constants shipped in the repository

It does not require a local GRAM installation or SPICE kernels.

Run the baseline no-GRAM examples:

```bash
julia --project=.AGORA examples/AGORA_Earth_NoGRAM.jl
julia --project=.AGORA examples/AGORA_Mars_NoGRAM.jl
```

The same onboarding mode is also available from the exported helper builders:

- `make_no_gram_planet`
- `make_no_gram_density_model`
- `make_no_gram_environment`

Minimal typed no-GRAM environment example:

```julia
using SpaceAGORA

env = make_no_gram_environment(
    planet=:earth,
    atmosphere=:none,
    EI_km=120.0,
    wind=false,
    topography=false,
)

env.ephemerides_model isa SimpleEphemeridesModel
# true
```

These presets are intended for onboarding, lightweight studies, and CI smoke coverage. They are not a replacement for GRAM- and SPICE-backed high-fidelity campaigns.

For full typed scenario construction and end-to-end `run_simulation(...)` examples, use the docs site rather than treating the README as a package manual.

## High-fidelity GRAM/SPICE mode

Use the GRAM/SPICE-backed path when you need:

- mission-quality atmospheric density and winds
- SPICE-based body orientation and geometry products
- third-body gravity
- solar-radiation-pressure ephemerides

High-fidelity runs expect a local GRAM installation rooted at `data/GRAMSuite.jl/GRAM Suite 2.0`. Treat `data/GRAMSuite.jl/GRAM Suite 2.0/Build/` as generated, host-specific output and rebuild it natively on each machine.

If you are part of the Space-FALCON Lab, use the lab-distributed GRAM assets. Otherwise, request GRAM from NASA here: [software.nasa.gov/software/MFS-33888-1](https://software.nasa.gov/software/MFS-33888-1).

For the current operational setup, keep the repository-root GRAM folder structure intact and ensure the required SPICE assets are available under the expected GRAM/SPICE tree.

## CLI/docs links

Package documentation:

- Docs site: [space-falcon-lab.github.io/SpaceAGORA.jl](https://space-falcon-lab.github.io/SpaceAGORA.jl)
- Local docs build: `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'`
- Build docs: `julia --project=docs docs/make.jl`

CLI entrypoint:

```bash
./bin/spaceagora
```

Common commands:

```bash
./bin/spaceagora run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1
./bin/spaceagora benchmark runtime-analysis smoke --output-dir=output/perf_cli
./bin/spaceagora assets check
```

For interactive scenario construction, asset policy, and the complete stable API surface, use the docs rather than internal modules under `src/`.

## Development notes

Generated runtime reports, telemetry CSVs, and local docs builds are not committed to git. They are regenerated into ignored directories such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.

Use the helper script below to regenerate the common local outputs:

```bash
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl docs
```

CI uploads telemetry CSVs as workflow artifacts instead of storing them in the repository.

`SpaceAGORA.run_simulation(...; isolate_state=true)` deep-copies the configuration by default so repeated or concurrent runs do not alias mutable state. Treat `isolate_state=false` as an expert-only performance lever.

An optional Docker/dev-container setup also exists for environment parity on some machines, but it is not required for first use and is not the primary onboarding path.
