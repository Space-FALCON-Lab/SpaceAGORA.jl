# Assets and Modes

Use this page when you need to decide whether your machine is ready to run
SpaceAGORA and which operating mode to choose.

This page is for users setting up local runs, deciding between no-GRAM and
higher-fidelity workflows, or checking whether required data is present.

Shortest successful command:

```text
julia --project=. src/cli/main.jl assets check
```

What to read next:

- [Quickstart](user/quickstart.md)
- [Installation & Environment](user/installation_environment.md)
- [First Simulation](user/first_simulation.md)

## Choose a mode

### No-GRAM baseline mode

Use this mode when you want onboarding, smoke coverage, or a lightweight local
run.

The default minimal quickstart path uses:

- `NoAtmosphereModel()`
- `SimpleEphemeridesModel()`

Optional open fallback atmosphere for no-GRAM runs:

- `ExponentialAtmosphereModel(planet)`
- `PiecewiseExponentialAtmosphereModel(h_breaks_m, rho_refs, Hs; ...)`
- `NRLMSISE00AtmosphereModel(; f107a, f107, ap, index_provider=..., use_space_indices=false)`

Notes:

- `ExponentialAtmosphereModel(planet)` is a single-scale-height, zero-wind approximation intended for a limited altitude band near its reference altitude. Its advisory validity range defaults to `[h_ref, h_ref + 5H]`.
- `PiecewiseExponentialAtmosphereModel(...)` is the multi-layer option for bounded studies that need a better altitude-shape fit while still avoiding GRAM assets.
- `NRLMSISE00AtmosphereModel(...)` accepts fixed geophysical indices, a callable provider, or `use_space_indices=true` to fetch CelesTrak F10.7 and Ap indices through `SpaceIndices`. Call `init_nrlmsise_space_indices!()` before long runs if you want the dataset prewarmed outside the solver. The documented standard NRLMSISE-00 validity band remains roughly `0 km` to `1000 km`.

It does not require:

- a GRAM installation
- local SPICE kernels
- GRAM surrogate grids

### High-fidelity mode

Use this mode when you need mission-quality atmosphere, SPICE-backed geometry,
or other higher-fidelity environment inputs.

It depends on machine-local or user-provided assets.

#### GRAM

Expected root:

```text
data/GRAMSuite.jl/GRAM Suite 2.0
```

This is a licensed external dependency. SpaceAGORA does not treat GRAM as part
of the baseline open onboarding path.

#### SPICE kernels

Expected under:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
```

These are used for high-fidelity frame orientation, SRP geometry, and N-body
ephemerides.

#### Gravity and topography harmonics

Expected under:

- `data/Gravity_harmonics_data`
- `data/Topography_harmonics_data`

These are optional and only required for studies or missions that enable those
models.

#### GRAM surrogate/static-grid bundle

Expected under:

```text
data/GRAM_surrogate
```

This bundle is optional. It accelerates selected GRAM-backed workflows but is
not required for the no-GRAM baseline mode.

## Check the current machine

The package-owned asset check distinguishes:

- baseline onboarding availability
- optional high-fidelity assets
- paths that are missing but non-blocking for no-GRAM mode

Useful commands:

```text
julia --project=. src/cli/main.jl assets check
julia --project=. src/cli/main.jl assets manifest
julia --project=. src/cli/main.jl assets setup-open
```

## Machine-readable manifest

Canonical manifest:

```text
data/assets_manifest.toml
```

This file is the machine-readable asset source of truth for:

- baseline built-in mode
- repository-owned open assets
- RPO station reference geometry, including the Gateway Core STL
- user-provided or licensed high-fidelity assets

You can inspect it with either:

```text
julia --project=. src/cli/main.jl assets manifest
julia --project=. scripts/assets/show_asset_manifest.jl
```

## Setup and check scripts

Baseline/open-mode bootstrap:

```text
julia --project=. scripts/assets/setup_open_assets.jl
```

This does not download GRAM or SPICE. It reports the baseline/open asset
contract and explicitly leaves licensed assets as user-provided.

Direct asset check script:

```text
julia --project=. scripts/assets/check_assets.jl
```
