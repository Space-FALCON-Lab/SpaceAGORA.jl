# Getting Started

## Clone and instantiate

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=.AGORA -e 'using Pkg; Pkg.instantiate()'
```

`/.AGORA` is the canonical committed execution environment for examples, benchmarks, tests, and CI. It is tracked in the repository, so there is no separate bootstrap step that copies the root `Project.toml` or `Manifest.toml` to create it at runtime.

## GRAM and data prerequisites

SpaceAGORA can run with local fallback atmosphere models, but high-fidelity GRAM-based studies require a local GRAM installation under `data/GRAMSuite.jl/GRAM Suite 2.0` and the expected SPICE assets.

For the current operational setup, use the repository root `data/GRAMSuite.jl/GRAM Suite 2.0` folder structure and keep host-specific GRAM build outputs local to the machine where GRAM is compiled.

## No-GRAM quickstart

SpaceAGORA now supports a first-class no-GRAM onboarding mode. This path uses:

- `NoAtmosphereModel()` or `ExponentialAtmosphereModel(planet)`
- `SimpleEphemeridesModel()`
- planet constants already shipped in the repository

It does not require a local GRAM installation or SPICE kernels.

Run the baseline no-GRAM examples directly:

```bash
julia --project=.AGORA examples/AGORA_Earth_NoGRAM.jl
julia --project=.AGORA examples/AGORA_Mars_NoGRAM.jl
```

The same mode is available from typed configuration builders through:

- `make_no_gram_planet`
- `make_no_gram_density_model`
- `make_no_gram_environment`

These presets are intended for onboarding, lightweight studies, and CI smoke coverage. They are not a replacement for GRAM- and SPICE-backed high-fidelity campaigns.

You can also exercise the stable root API without GRAM assets:

```jldoctest; setup = :(using SpaceAGORA)
julia> profile = parse_parallel_profile("R2")
R2

julia> cfg = simulation_engine_config_from_env(Dict(
           "SPACEAGORA_PARALLEL_PROFILE" => parallel_profile_name(profile),
           "SPACEAGORA_SAVE_BUNDLE" => "0",
       ));

julia> (cfg.parallel.profile, cfg.artifacts.save_bundle)
("R2", false)
```

## High-fidelity GRAM mode

Use the GRAM and SPICE-backed path when you need:

- mission-quality atmospheric density and winds
- third-body gravity
- solar-radiation pressure ephemerides
- SPICE-based body orientation and geometry products

That path still expects a local GRAM installation under `data/GRAMSuite.jl/GRAM Suite 2.0` together with the required SPICE assets.

## Run a minimal example

```bash
julia --project=.AGORA examples/Earth_Thruster_Test.jl
```

## Run the telemetry verification study

```bash
julia --project=.AGORA benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

## Runtime performance knobs

`SpaceAGORA.run_simulation(...; isolate_state=true)` deep-copies the configuration by default so repeated or concurrent runs cannot alias mutable state. For large mission definitions or throughput-heavy short-run sweeps, expert callers can evaluate `isolate_state=false` as a setup-cost lever, but only when they own the configuration instance and will not reuse it across overlapping or state-mutating runs.

The runtime analysis benchmark reports this copy/setup overhead explicitly, so use those results before turning isolation off for a workflow.

## Generated outputs stay local

Generated runtime reports, telemetry CSVs, and docs builds are regenerated into ignored directories such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`. This keeps machine-local absolute paths, usernames, and host-specific build state out of the repository.

Use the helper script when you want to recreate the common local artifacts:

```bash
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl docs
```

## Build the documentation locally

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
