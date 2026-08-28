# SpaceAGORA.jl

SpaceAGORA.jl is a Julia toolkit for spacecraft orbit and atmospheric-flight
simulation, including aerobraking, entry, and aerocapture workflows. The main
user-facing surfaces are:

- the root `SpaceAGORA` Julia API
- the repository CLI wrapper at `./bin/spaceagora` on Linux/macOS or
  `bin\spaceagora.bat` on Windows
- runnable example scripts under `examples/`

Cross-platform command conventions used below:

- `julia --project=. src/cli/main.jl ...` works on Windows, Linux, and macOS
- `julia --project=. scripts/ensure_gram_native.jl ...` works on Windows,
  Linux, and macOS
- `./bin/spaceagora` and `./scripts/ensure_gram_native.sh` remain available as
  convenience wrappers on Linux/macOS

## Installation

Use the repository root environment as the canonical committed execution environment for examples, tests, and normal local runs:

```text
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

That gives you the baseline open-data onboarding path immediately. You do not
need GRAM or SPICE to run the first quickstart example, and there is no bootstrap step beyond instantiating the committed root environment.

Helpful first checks:

```text
julia --project=. src/cli/main.jl assets check
julia --project=. src/cli/main.jl assets setup-open
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
step-by-step local setup guide for licensed GRAM assets. The full GRAM submodule
checkout uses Git LFS and downloads several GB of binary data, so it is separate
from the default baseline installation path.

### Building the native GRAM library

GRAM is a C/C++ codebase that SpaceAGORA calls through FFI, so the shared
library has to be compiled once per machine. One command does all of it:

```text
julia --project=. scripts/ensure_gram_native.jl
```

That returns immediately if the library for your host is already there, so it
is safe to run unconditionally before a GRAM-backed run. Otherwise it builds
`data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib/libGRAM.so` (`.dylib` on macOS,
`.dll` on Windows), which takes on the order of 20 seconds on 24 cores.

You need a C/C++ toolchain and GNU Make (`make` on Linux, `gmake` on macOS via
`brew install make`, `mingw32-make` from MSYS2 on Windows). You do **not** need
a Fortran compiler, and you do not need to install CSPICE on Linux x86_64 or
Windows — a suitable archive is bundled and selected automatically. On macOS,
run `brew install cspice` first.

Force a full rebuild with:

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

You need `--clean` in one specific situation: **the GRAM tree was copied or
`rsync`ed in from another machine with its `Build/lib` already populated** — as
happens with `scripts/remote/spaceagora-remote`, which mirrors the working
tree. The command above skips the build whenever `Build/lib/libGRAM.<ext>`
exists, and a foreign library satisfies that check, so you end up loading a
binary built for the wrong host. Moving or renaming the checkout on the *same*
machine is detected and handled automatically; no flag needed.

#### Build prerequisites

| Requirement | Notes |
|---|---|
| C/C++ toolchain | GCC on Linux; Clang via Xcode command-line tools on macOS; MinGW-w64 (MSYS2) on Windows |
| GNU Make | Linux: `make`. macOS: `gmake`, from `brew install make` — the GRAM makefiles use GNU Make features the `make` shipped with macOS does not support. Windows: `mingw32-make`. |
| CSPICE | Bundled and selected automatically on Linux x86_64 and Windows MinGW. macOS: `brew install cspice`. Linux arm64: supply your own archive. |
| Git LFS | The GRAM submodule stores its surrogate payloads and SPICE kernels via LFS. |

No Fortran compiler is required, despite `gfortran` appearing in GRAM's build
configuration — `Build/makefile.defs` ends with `undefine FC`, so the library
build never invokes it.

#### What the command does

`ensure_gram_native.jl` delegates to the vendored helper
`GRAM Suite 2.0/simulation/GRAM/build_gram.sh`, which runs
`Build/setup_cspice.sh` (you never run this yourself), then `make shared` with
one job per core, then writes `simulation/GRAM/gram.env` and
`simulation/GRAM/.gram-build-manifest` recording the host and root path the
artifacts belong to.

`shared` is the only make target that produces `libGRAM`. Plain `make` and the
per-planet targets such as `make Mars` build static libraries and leave
`Build/lib` without the shared library SpaceAGORA loads, so there is no faster
single-planet route.

`gram.env` and `.gram-build-manifest` hold absolute paths for the machine that
wrote them. Never commit them — a committed one makes the next person's build
fail against a directory that does not exist on their machine.

#### If the build fails

| Symptom | Cause and fix |
|---|---|
| `Build finished but shared library not found`, naming a path that is not yours | A `gram.env`/`.gram-build-manifest` from another machine is present. Delete both and re-run. |
| `Could not auto-find a Linux CSPICE archive` | Linux arm64 has no bundled archive. Install CSPICE so `/usr/lib/libcspice.a` (or `/usr/local/lib`, `/usr/lib64`, `/usr/lib/aarch64-linux-gnu`) exists, or build with `make shared -j SPICE_LIB=/absolute/path/to/cspice.a` from `data/GRAMSuite.jl/GRAM Suite 2.0/Build`. |
| `gmake not found` on macOS | `brew install make`. |
| `make` succeeded but there is still no `libGRAM` | You ran `make` or a per-planet target instead of `make shared`. Prefer `ensure_gram_native.jl`, which always builds the right target. |
| `GRAMAtmosphereModel` load error although the library exists | It was built on a different machine — rebuild with `--clean`. |

The full walkthrough, including copying the licensed NASA GRAM folders into
place, is in
[docs/src/user/gramsuite_setup.md](docs/src/user/gramsuite_setup.md).

## Quick Start

The fastest first run is:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

Or through the CLI:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl
```

That example is intentionally simple:

- Earth orbit
- no atmosphere drag
- no SPICE dependency
- no GRAM dependency
- results written under `output/`

If you only want to verify the end-to-end path with a short run, use:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --smoke
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

```text
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
  `julia --project=. scripts/ensure_gram_native.jl`

Useful commands:

```text
julia --project=. src/cli/main.jl assets check
julia --project=. scripts/ensure_gram_native.jl
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

If you want fuller mission examples after that, start with:

- `examples/AGORA_Vex.jl`
- `examples/AGORA_Odyssey.jl`

Both of those use the same GRAM access pattern as the basic GRAM example, but
add mission-specific guidance, maneuvers, and higher-fidelity dynamics.

## CLI and docs

- Docs site: [space-falcon-lab.github.io/SpaceAGORA.jl](https://space-falcon-lab.github.io/SpaceAGORA.jl)
- Cross-platform CLI entrypoint: `julia --project=. src/cli/main.jl`
- Linux/macOS wrapper: `./bin/spaceagora`
- Windows wrapper: `bin\spaceagora.bat`
- Run an example: `julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl`
- Asset check: `julia --project=. src/cli/main.jl assets check`
- Verification study: `julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1`

## For contributors

- User docs live in `docs/src/` and are built with `julia --project=. docs/make.jl`
- Generated API pages come from `docs/public_api_symbols.jl`
- Architecture and quality references live under the docs Maintainer Guide

Generated reports, docs builds, and other local outputs stay in ignored paths
such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.
