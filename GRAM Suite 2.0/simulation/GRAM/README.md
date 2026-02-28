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

```bash
julia "GRAM Suite 2.0/simulation/GRAM/build_offline_static_grids.jl" \
  --planets=all \
  --out-dir="GRAM Suite 2.0/simulation/GRAM/static_grids/p175_adaptive_all_planets" \
  --surrogate-only=true \
  --surrogate-source=direct \
  --surrogate-lat-step-deg=1.75 \
  --surrogate-lon-step-deg=1.75 \
  --alt-min-km=0 \
  --alt-max-km=300 \
  --surrogate-adaptive-alt=true \
  --surrogate-adaptive-min-step-km=0.5 \
  --surrogate-adaptive-max-step-km=6.0 \
  --surrogate-adaptive-pilot-alt-step-km=2.0 \
  --surrogate-adaptive-pilot-lat-step-deg=10.0 \
  --surrogate-adaptive-pilot-lon-step-deg=10.0
```
