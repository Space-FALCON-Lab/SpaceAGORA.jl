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
