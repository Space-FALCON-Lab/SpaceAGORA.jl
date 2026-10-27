# Assets and Modes

A simulation's environment inputs decide which assets it needs. The no-GRAM
baseline, used for first runs, needs nothing beyond the repository. Native
GRAM, SPICE-backed geometry and the surrogate workflows need the additional
inputs described below.

After [preparing the core environment](user/installation_environment.md), list
what this machine has by running, from the repository root:

```text
julia --project=. src/interfaces/cli/main.jl assets check
```

The check also covers the reference data the verification studies use. A
missing required reference makes the check fail even though the no-GRAM
[Quickstart](user/quickstart.md) still runs, so read the asset names in the
report against the workflow you have in mind. [Your First Script](user/first_script.md)
explains the baseline configuration.

## Choose a mode

### No-GRAM baseline mode

This mode is for onboarding, smoke runs and lightweight local runs. The
quickstart uses:

- `NoAtmosphereModel()`
- `SimpleEphemeridesModel()`

It needs no GRAM installation, no local SPICE kernels and no GRAM surrogate
grids. When a no-GRAM run needs an atmosphere, three open models are available:

- `ExponentialAtmosphereModel(planet)` is a single-scale-height model with no
  wind, meant for a limited altitude band near its reference altitude. Its
  advisory validity range defaults to `[h_ref, h_ref + 5H]`.
- `PiecewiseExponentialAtmosphereModel(h_breaks_m, rho_refs, Hs; ...)` is the
  multi-layer version, for bounded studies that need a better fit to the
  altitude profile without GRAM assets.
- `NRLMSISE00AtmosphereModel(; f107a, f107, ap, index_provider=..., use_space_indices=false)`
  takes fixed geophysical indices, a callable provider, or
  `use_space_indices=true` to fetch CelesTrak F10.7 and Ap indices through
  `SpaceIndices`. To load that dataset before a long run rather than during the
  solve, call `init_nrlmsise_space_indices!()` first. The model's advisory
  validity range is 0 to 1000 km, the standard NRLMSISE-00 band.

### High-fidelity mode

This mode is for models that need native GRAM atmosphere data, SPICE-backed
geometry or other inputs that live on the local machine or that you provide.

| Asset | Expected location | Notes |
|---|---|---|
| GRAM | `data/GRAMSuite.jl/GRAM Suite 2.0` | a licensed external dependency, not part of the baseline onboarding path |
| SPICE kernels | `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE` | frame orientation, SRP geometry and n-body ephemerides |
| Gravity and topography harmonics | `data/planets/gravity`, `data/planets/topography` | needed only by studies or missions that enable those models |
| GRAM surrogate grids | `data/planets/gram_surrogate` | optional; used by selected GRAM-backed workflows, never by the no-GRAM baseline |

### What a fresh clone contains

After `GIT_LFS_SKIP_SMUDGE=1 git clone --filter=blob:none`, the checkout holds
the following.

| Asset | In the clone? | Access | Needed by |
|---|---|---|---|
| Built-in planets, fallback atmospheres, `SimpleEphemeridesModel` | yes | open | the quickstart, `AGORA_Earth_NoGRAM.jl`, `AGORA_Earth_MonteCarlo.jl` |
| Gravity harmonics CSVs (`data/planets/gravity`) | yes | open | missions and the verification study |
| Topography harmonics (`data/planets/topography`) | yes | open | nothing on the default path |
| Odyssey, VEx and GMAT Earth truth files (`data/telemetry/...`, under 3 MB) | yes | open | the verification study |
| RPO Gateway Core STL (`data/scene/models/gateway_core.stl`) | yes | open | the RPO examples |
| GRAM surrogate grids (`data/planets/gram_surrogate/*.jls`, seven Git LFS objects, about 2.6 GB) | pointers only with `GIT_LFS_SKIP_SMUDGE=1`, downloaded otherwise | Git LFS | one benchmark study; at run time, surrogates are resolved through the submodule |
| `data/GRAMSuite.jl` submodule: the vendored `GRAMSuite` package, the SPICE kernels and the GRAM Suite 2.0 tree | empty directory until you initialize it | public Git repository with Git LFS payloads of several GB; [GRAMSuite Setup](user/gramsuite_setup.md) says what is fetched and what you must provide | GRAM-backed examples, `Earth("", SPICE_PATH)`-style planets, the verification study and the benchmark launchers |
| Simulator reference sets (Basilisk, GMAT parity matrix) and public OPS-SAT-1 telemetry | no (Git ignores them) | lab organization, `scripts/dev/fetch_private_telemetry.sh references` | the parity matrix `test/gmat_scenario_matrix.jl` |
| CYGNSS and GRIFEX flight telemetry | no (Git ignores them) | per person, granted by the PI; never commit or share it further | the mission studies that name them |

`assets check` reports the same picture for the machine you are on. Each entry
shows one state and its licensing class:

| State | Meaning |
|---|---|
| `usable` | a file holds real (non-pointer) bytes, or a directory has at least one nonempty payload and no pointers; built-in entries are always usable |
| `missing-required`, `missing-optional` | the expected file or directory does not exist |
| `empty-required`, `empty-optional` | the directory has no visible files |
| `pointer-only-required`, `pointer-only-optional` | some or all files are undownloaded Git LFS pointers |
| `invalid-required`, `invalid-optional` | the file is empty (0 bytes), or every visible file in the directory is |

Every entry that is not usable gets a `reason` line, and the report ends with a
summary. The command exits with status 1 when a required entry is not usable
and 0 otherwise; optional entries only affect the report. An error while
inspecting the file system, such as a missing read permission, stops the check
instead.

For example, after a clone with `GIT_LFS_SKIP_SMUDGE=1`,
`gram_surrogate_directory` reports `pointer-only-optional` while the seven
`.jls` files are still pointers; `git lfs ls-files` marks them with `-` instead
of `*`, and `git lfs pull --include "data/planets/gram_surrogate/*"` fetches
them before a surrogate-backed benchmark. The reason line names the Git
checkout that owns each pointer, a GRAMSuite submodule or linked asset
included. On Unix it gives a quoted shell command; on Windows, where quoting
depends on the shell, it gives the owning checkout and the Git LFS include
pattern instead.

The check skips hidden files and directories. It follows file symlinks, and an
asset root that is itself a directory symlink, but does not descend into
nested directory symlinks. It reads directory listings, sizes and small file
headers only, and never downloads anything. `usable` therefore means bytes are
present, not that the data is complete, scientifically valid or ready for
native use; a directory that holds both empty and nonempty files can be
`usable`. To validate a grid's contents, use `assets check --preset` (see
[Named frozen surrogate presets](#Named-frozen-surrogate-presets)).

### Verification reference data

The simulator-to-simulator reference sets that the parity matrix
(`test/gmat_scenario_matrix.jl`) grades against, and the flight telemetry the
mission studies use, are not tracked in this repository.
`scripts/dev/fetch_private_telemetry.sh` copies them into ignored
`data/telemetry/` paths; `data/telemetry/PRIVATE_TELEMETRY.md` gives the
sources, the access rules and what happens without them (the affected test
sets are skipped).

## Asset commands and the manifest

`data/assets_manifest.toml` is the machine-readable list of assets: the
built-in baseline, the open assets in the repository, the RPO station reference
geometry including the Gateway Core STL, and the user-provided or licensed
high-fidelity assets. Each entry is marked required or optional. The telemetry
truth files are required, so without them the asset check fails, although a
no-GRAM simulation still runs.

```text
julia --project=. src/interfaces/cli/main.jl assets check
julia --project=. src/interfaces/cli/main.jl assets manifest
julia --project=. src/interfaces/cli/main.jl assets setup-open
```

`assets manifest` prints the manifest. `assets setup-open` lists the open
entries and states that the licensed ones remain yours to provide; it downloads
neither GRAM nor SPICE. The same three operations exist as scripts:

```text
julia --project=. scripts/assets/check_assets.jl
julia --project=. scripts/assets/show_asset_manifest.jl
julia --project=. scripts/assets/setup_open_assets.jl
```

## Named frozen surrogate presets

For the Odyssey guidance and control exercise, follow the
[surrogate tutorial](tutorials/odyssey_surrogate.md). Its environment,
`examples/odyssey_surrogate_env`, retrieves a pinned public GRAMSuite wrapper;
it needs neither native GRAM nor a private repository. The published preset
grids and payloads are CC BY 4.0: cite the preset by name and version, and
credit NASA's GRAM Suite (Mars-GRAM). The near-surface payload also holds
Mars-GRAM's MOLA terrain values at its lattice nodes, so credit NASA MOLA as its
archive README asks. The scenario assets keep their original NAIF and
gravity-data terms.

```sh
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets list
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets fetch --preset=odyssey_p20_frozen_v1 --version=1.0.0
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets check --preset=odyssey_p20_frozen_v1 --version=1.0.0
```

`fetch` installs the payload and verifies its exact bytes. `check` stays offline
and also validates the schema, the coordinates and the generation metadata. The
example retrieves its preset the same way when it starts. A file given with
`--file` must match the named preset and never falls back to a download.
Missing files, LFS pointers, corrupt downloads and incompatible data are errors
before any propagation starts. Julia's artifact system handles concurrent
installation; a corrupt managed artifact or override must be repaired by hand,
and a file you supplied is never replaced.

A named version fixes the atmosphere's epoch and domain, and changing the
simulation date does not move it. The catalog lists four presets, all frozen at
the P20 instant (see [Atmosphere Models](user/atmosphere_models.md)):

- `odyssey_p20_frozen_v1`, the bounded Odyssey grid;
- `mars_global_upper_p20_frozen_v1`, a global Mars grid from 80 to 365 km;
- `mars_global_near_surface_p20_frozen_v1`, a near-surface Mars atmosphere from
  5 m above the surface. Versions 1.0.0 and 1.1.0 store no winds; 1.2.0 adds
  them. Version 1.0.0 reaches 75 km, and 1.1.0 and 1.2.0 reach 81 km areoid
  height, above the global grid's 80 km floor everywhere;
- `mars_global_descent_p20_frozen_v1`, a composite with no data of its own: the
  near-surface preset 1.2.0 below 80 km above the ellipsoid and the global
  preset 1.0.0 at and above it, for descents from orbit to the ground.

The same `fetch` and `check` commands take
`--preset=mars_global_upper_p20_frozen_v1 --version=1.0.0`, whose archive is
about 230 MB, or `--preset=mars_global_near_surface_p20_frozen_v1 --version=1.2.0`
(or `--version=1.1.0` or `--version=1.0.0`), whose archives are about 34 MB,
17 MB and 15 MB. For `--preset=mars_global_descent_p20_frozen_v1 --version=1.0.0`,
`fetch` installs and verifies both components (the upper 1.0.0 and near-surface
1.2.0 archives) and `check` builds the composite from them; it takes no
`--file`. Other planets and epochs need separately generated and validated
presets. The results bundle's manifest has an `atmosphere` section that
identifies the backend, the preset and version, the payload and catalog hashes
and the frozen-domain contract.
