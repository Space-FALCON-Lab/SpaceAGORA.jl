# GRAM Simulation Drop-In

This folder is a portable bridge for using GRAM from a simulation project.

## What this gives you

- `build_gram.sh`: Bash build script for macOS/Linux/MSYS shells.
- `run_julia_smoke_test.sh`: Bash Julia smoke test.
- `build_gram.ps1`: native Windows PowerShell build script.
- `run_julia_smoke_test.ps1`: native Windows PowerShell smoke test.
- `build_gram.cmd` and `run_julia_smoke_test.cmd`: Windows CMD wrappers.
- `julia_smoke_test.jl`: the actual Julia smoke test.
- `gram.env` / `gram.env.ps1`: generated env files with `GRAM_ROOT` and `GRAM_LIB`.

## Quick start

macOS/Linux/MSYS shell:

1. `./build_gram.sh /absolute/path/to/GRAM Suite 2.0`
2. `./run_julia_smoke_test.sh`

Windows PowerShell:

1. `.\build_gram.ps1 "C:\path\to\GRAM Suite 2.0"`
2. `.\run_julia_smoke_test.ps1`

Windows CMD:

1. `build_gram.cmd "C:\path\to\GRAM Suite 2.0"`
2. `run_julia_smoke_test.cmd`

## Notes

- macOS uses `gmake` (GNU make). Linux uses `make`.
- Windows uses `mingw32-make` (or `make`) from MSYS2/MinGW.
- Scripts call `Build/setup_cspice.sh` automatically.
- On Windows, bundled `common/cspice/lib/cspice_mingw64.a` is used.

## Offline Atmosphere Products

Use `build_offline_static_grids.jl` to generate offline GRAM payloads (`*.jls`) for SpaceAGORA.

Default surrogate profile (`p175_mid`, recommended):

```bash
julia "GRAM Suite 2.0/simulation/GRAM/build_offline_static_grids.jl" \
  --planets=all \
  --out-dir="GRAM Suite 2.0/simulation/GRAM/static_grids/p175_mid_all_planets" \
  --surrogate-only=true \
  --surrogate-source=direct \
  --surrogate-lat-step-deg=1.75 \
  --surrogate-lon-step-deg=1.75 \
  --alt-min-km=0 \
  --alt-max-km=300
```

Earth low-altitude caveat:

- Earth GRAM can assert in some MERRA2 conditions below ~5 km.
- If that occurs, build Earth separately at `5-300 km`:

```bash
julia "GRAM Suite 2.0/simulation/GRAM/build_offline_static_grids.jl" \
  --planets=earth \
  --out-dir="GRAM Suite 2.0/simulation/GRAM/static_grids/p175_mid_all_planets" \
  --surrogate-only=true \
  --surrogate-source=direct \
  --surrogate-lat-step-deg=1.75 \
  --surrogate-lon-step-deg=1.75 \
  --alt-min-km=5 \
  --alt-max-km=300 \
  --month=1 \
  --earth-merra2-path="GRAM Suite 2.0/Earth/data/MERRA2data"
```

Adaptive altitude grid (dynamic precision):

- Enable `--surrogate-adaptive-alt=true` to automatically densify altitude nodes where density/gradient is high and coarsen where atmosphere is smoother.
- Keep latitude/longitude spacing fixed; altitude spacing becomes nonuniform per-planet.
- Frozen per-planet policy is enabled by default (`--use-frozen-planet-policy=true`).

Frozen policy (`2026-02-28_rho1e-13_relp95-0p10_vmax250_ceiling-m365-v460-e1115`):

- Earth: `alt=5..1115 km`, `min_step=0.5 km`, `max_step=250 km`
- Mars: `alt=0..365 km`, `min_step=0.5 km`, `max_step=50 km`
- Venus: `alt=0..460 km`, `min_step=0.5 km`, `max_step=30 km`
- Jupiter: `alt=0..1000 km`, `min_step=0.5 km`, `max_step=250 km`
- Titan: `alt=0..2500 km`, `min_step=0.5 km`, `max_step=200 km`
- Neptune: `alt=0..4000 km`, `min_step=0.5 km`, `max_step=250 km`
- Uranus: `alt=0..7000 km`, `min_step=0.5 km`, `max_step=250 km`

Earth month default:

- Builder now defaults Earth to `month=1` when `--month` is not explicitly set (to avoid known MERRA/GRAM dew-point assertion cases seen for some months during dense global sweeps).
- To force a specific month for all planets, pass `--month=<1..12>`.

```bash
julia "GRAM Suite 2.0/simulation/GRAM/build_offline_static_grids.jl" \
  --planets=all \
  --out-dir="GRAM Suite 2.0/simulation/GRAM/static_grids/p175_adaptive_all_planets" \
  --surrogate-only=true \
  --surrogate-source=direct \
  --surrogate-lat-step-deg=1.75 \
  --surrogate-lon-step-deg=1.75 \
  --surrogate-adaptive-alt=true \
  --surrogate-adaptive-pilot-alt-step-km=2.0 \
  --surrogate-adaptive-pilot-lat-step-deg=10.0 \
  --surrogate-adaptive-pilot-lon-step-deg=10.0
```

Overrides:

- Disable frozen policy: `--use-frozen-planet-policy=false`
- Global altitude range: `--alt-min-km`, `--alt-max-km`
- Per-planet altitude range: `--alt-min-km-<planet>`, `--alt-max-km-<planet>`
- Global adaptive steps: `--surrogate-adaptive-min-step-km`, `--surrogate-adaptive-max-step-km`
- Per-planet adaptive steps: `--surrogate-adaptive-min-step-km-<planet>`, `--surrogate-adaptive-max-step-km-<planet>`

## Runtime policy (SpaceAGORA)

- Default runtime policy is now: **offline surrogate first**, then **fallback to point-to-point GRAM** when outside surrogate coverage or when unsupported feature knobs are enabled.
- Default surrogate location:
  - `GRAM Suite 2.0/<Planet>/<planet>_surrogate.jls` (for example: `GRAM Suite 2.0/Mars/mars_surrogate.jls`)
  - Legacy fallback is still supported: `GRAM Suite 2.0/simulation/GRAM/static_grids/p175_mid_all_planets/surrogates`

Environment knobs:

- `SPACEAGORA_GRAM_OFFLINE_SURROGATE=on|off|auto` (default: `on`)
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_DIR=/abs/path/to/surrogates`
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_DIR_<PLANET>=/abs/path/to/surrogates`
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_FILE=/abs/path/to/<planet>_surrogate.jls`
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_FILE_<PLANET>=...` (e.g., `..._MARS`)
- `SPACEAGORA_GRAM_WIND_MODE=auto|nominal|perturbed` (default: `auto`)
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M=<meters>`
- `SPACEAGORA_GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M_<PLANET>=<meters>`

Notes:

- `SPACEAGORA_GRAM_OFFLINE_SURROGATE=auto` means:
  - behave like `on` if surrogate use is supported for the active GRAM feature configuration
  - otherwise fall back to point-to-point GRAM
- Above surrogate altitude ceiling, runtime returns `rho=0` (vacuum behavior) for that sample.
- Titan uses a default hybrid policy for stability: surrogate above `260 km`, point-to-point GRAM below `260 km`.

Legacy static-grid mode remains available, but the default path uses offline surrogate payloads.
