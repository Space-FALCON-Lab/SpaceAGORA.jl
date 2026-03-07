# Data and Assets

SpaceAGORA now supports two operational modes:

## 1. No-GRAM baseline mode

This path is intended for onboarding and CI smoke coverage.

It uses:

- `NoAtmosphereModel()`
- `ExponentialAtmosphereModel(planet)`
- `SimpleEphemeridesModel()`

It does not require:

- a GRAM installation
- local SPICE kernels
- GRAM surrogate grids

## 2. High-fidelity mode

This path depends on machine-local or user-provided assets.

### GRAM

Expected root:

```text
data/GRAMSuite.jl/GRAM Suite 2.0
```

This is a licensed external dependency. SpaceAGORA does not treat GRAM as part of the baseline open onboarding path.

### SPICE kernels

Expected under:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
```

These are used for high-fidelity frame orientation, SRP geometry, and N-body ephemerides.

### Gravity and topography harmonics

Expected under:

- `data/Gravity_harmonics_data`
- `data/Topography_harmonics_data`

These are optional and only required for studies or missions that enable those models.

### GRAM surrogate/static-grid bundle

Expected under:

```text
data/GRAM_surrogate
```

This bundle is optional. It accelerates selected GRAM-backed workflows but is not required for the no-GRAM baseline mode.

## Asset check

Use the package-owned CLI to inspect the current machine state:

```bash
./bin/spaceagora assets check
```

That command distinguishes:

- baseline onboarding availability
- optional high-fidelity assets
- paths that are missing but non-blocking for no-GRAM mode

## Machine-readable manifest

Canonical manifest:

```text
data/assets_manifest.toml
```

This file is the machine-readable asset contract for:

- baseline built-in mode
- repository-owned open assets
- user-provided or licensed high-fidelity assets

You can inspect it with either:

```bash
./bin/spaceagora assets manifest
julia --project=.AGORA scripts/assets/show_asset_manifest.jl
```

## Setup and check scripts

Baseline/open-mode bootstrap:

```bash
julia --project=.AGORA scripts/assets/setup_open_assets.jl
```

This does not download GRAM or SPICE. It reports the baseline/open asset contract and
explicitly leaves licensed assets as user-provided.

Direct asset check script:

```bash
julia --project=.AGORA scripts/assets/check_assets.jl
```
