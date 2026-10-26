# GRAMSuite Setup

This page is for users who have been granted access to NASA's GRAM Suite and
want to run the native GRAM-backed models in a local SpaceAGORA checkout. The
baseline no-GRAM path does not need any of it, and neither do the grid
surrogates described on [GRAM Grid Surrogates](gram_grid_surrogates.md).

Related pages:

- [Installation & Environment](installation_environment.md)
- [Assets & Modes](../assets.md)
- [Verification Study](verification_study.md)

## What this document covers

Native GRAM in SpaceAGORA comes from two places. The public `GRAMSuite.jl`
wrapper is a Git submodule of this repository at `data/GRAMSuite.jl`. It holds
the Julia integration and, in a folder named `GRAM Suite 2.0`, the SPICE
kernels and the native build scripts. The GRAM Suite itself, the source code
and data that the native library is built from, is NASA software that you
request from NASA. It is not in the public checkout, and this page does not
distribute it.

The setup has four steps: pull the submodule, request and receive the NASA
distribution, copy it into the submodule's `GRAM Suite 2.0` folder, and check
that SpaceAGORA finds the GRAM and SPICE directories. Building the native
library and a first GRAM-backed run complete it.

## Request GRAM access

Request access through the NASA Software Catalog entry for GRAM:

- <https://software.nasa.gov/software/MFS-33888-1>

Once NASA approves the request, obtain the GRAM Suite package the way NASA
specifies. Keep the original delivery intact until you have finished copying
and checking it as described below.

## Pull the `GRAMSuite.jl` submodule

Work from the SpaceAGORA repository root. The wrapper stores its SPICE kernels
and static surrogate grids in Git LFS, several gigabytes in all, so before you
start make sure the machine has room for the final files and for Git LFS's
temporary downloads. In practice that means 15 to 20 GB of free space.

If you already have a wrapper checkout with local work in it, keep that
checkout and use a separate one for these steps; do not reset it.

```text
git submodule update --init --recursive -- data/GRAMSuite.jl
```

This checks out the exact submodule commit that SpaceAGORA records. Do not add
`--remote`, which selects a branch tip instead of that commit.

Check the recorded commit against the wrapper's actual commit, and look for
local changes:

```text
git rev-parse HEAD:data/GRAMSuite.jl
git -C data/GRAMSuite.jl rev-parse HEAD
git -C data/GRAMSuite.jl status --short
```

For the public wrapper the two commit identifiers must agree. The preparation
helper in the next section also checks the source bytes and the public pin in
`.github/gramsuite-revisions`; a matching folder name is not enough.

After the checkout the wrapper is at:

```text
data/GRAMSuite.jl
```

The wrapper is not the GRAM Suite. It gives SpaceAGORA the Julia integration
layer and the expected folder layout; the NASA distribution supplies the GRAM
content that goes into that layout.

To look at the wrapper source or the layout without downloading the LFS files
yet, skip LFS smudging during the checkout:

```text
GIT_LFS_SKIP_SMUDGE=1 git submodule update --init --recursive -- data/GRAMSuite.jl
```

Later, once there is enough free disk space, fetch the LFS files from inside
the submodule:

```text
cd data/GRAMSuite.jl
git lfs pull
```

## Prepare package-loading dependencies offline

This optional step checks, without network access, that the packages load.
It applies when the pinned wrapper source and SpaceAGORA's package depot are
already on the machine. A wrapper checked out as source only can lack its own
dependency metadata even when SpaceAGORA itself imports, and that matters when
a loading helper, such as the examples' `setup_gram_example!`, puts the
wrapper at the front of Julia's load path.

From the SpaceAGORA repository root, with Python 3.11 or later and Git
available:

```sh
python3.11 scripts/prepare_gram_loading_environment.py \
  --gram-source "$PWD/data/GRAMSuite.jl" \
  --output "$HOME/spaceagora-gram-loading-prepared"
```

The output directory must be new, its parent must exist, and it must lie
outside both checkouts. The helper checks `Project.toml` and every file under
`src/` against the configured GRAM commit, using local Git objects. For a
source-only copy kept elsewhere, pass its path as `--gram-source` and add
`--gram-object-repo` pointing to a local GRAM repository that contains that
commit. The helper never fetches, so missing objects are an unmet
prerequisite. `--mode dev` selects the separately configured development pin,
which needs its own verified source and local objects. A result for the public
wrapper says nothing about the development wrapper.

The prepared directory holds the wrapper's unchanged `Project.toml`, a link to
its source, a dependency manifest, and `preparation.json` with the source and
metadata hashes. The manifest copies complete stanzas, with exact versions,
UUIDs and tree hashes, from SpaceAGORA's committed root manifest. It includes
the required transitive dependencies and list-form weak dependencies, whose
stanzas Julia needs to identify extension triggers. Dictionary-form weak
dependencies already carry their UUIDs and are left as they are. The root's
`project_hash` is left out because it describes SpaceAGORA's project, not the
wrapper's.

The helper refuses locally changed root project, manifest or pin files,
modified or extra wrapper source files, missing or ambiguous dependency
entries, and path dependencies. It does not resolve, instantiate or install
packages, import Julia packages, copy licensed assets, or modify either
checkout. Because the prepared directory links to the original source, keep
that source in place and unchanged while you use the preparation, and keep
`preparation.json` as its record.

Then run the separate loading diagnostic against the prepared directory:

```sh
JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_PKG_OFFLINE=true \
  julia --startup-file=no --compiled-modules=existing \
  test/smoke/gram_package_readiness.jl \
  --gram-project "$HOME/spaceagora-gram-loading-prepared" \
  --output "$HOME/spaceagora-gram-loading-report"
```

The report directory must also be new. The preparation supplies metadata, not
package files, so dependencies that are not installed still make the
diagnostic fail.
[GRAM package loading readiness](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/quality/gram_package_readiness.md)
describes the six import and failure cases and the runtime limits. A passing
report covers the package and extension identities it observed. It says
nothing about production helpers, workers, native GRAM libraries or data, or
simulation results. The native setup below is a separate step.

## Expected target location

SpaceAGORA expects the GRAM Suite tree at:

```text
data/GRAMSuite.jl/GRAM Suite 2.0
```

`scripts/ensure_gram_native.jl`, the examples' `SPICE_PATH` and the asset
report all use this path. The parts SpaceAGORA uses most are:

- `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE`
- `data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib`
- the rest of the GRAM Suite runtime tree that `GRAMSuite.jl` needs

The wrapper checkout already provides `SPICE`, `simulation` and `Julia` in
that folder. `Build` and the planet folders come from the NASA delivery.

## Copy the official GRAM Suite folders

Once NASA has given you the GRAM Suite files, copy the contents of the
distribution into that folder:

1. Unpack or open the GRAM Suite delivery from NASA.
2. Find its top-level `GRAM Suite 2.0` directory.
3. Copy the contents of that directory into `data/GRAMSuite.jl/GRAM Suite 2.0`.
4. Check that the final path is exactly `data/GRAMSuite.jl/GRAM Suite 2.0`.

After the copy, `Build` sits next to the folders the wrapper provides:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/Build
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
data/GRAMSuite.jl/GRAM Suite 2.0/simulation
```

Do not rename the `GRAM Suite 2.0` directory; SpaceAGORA and the wrapper
scripts expect exactly that name.

## Build or verify the native GRAM library

SpaceAGORA calls GRAM through a native shared library, `libGRAM`. Build it on
the machine that will use it: a library copied from another machine may not
work on this host, even when its filename has the expected extension.

### Build prerequisites

Use the compiler and CSPICE archive that suit your GRAM distribution and host
architecture; the wrapper's `data/GRAMSuite.jl/README.md` lists them for each
platform. The Unix build helper uses GNU Make: `make` on Linux and `gmake` on
macOS, where `brew install make` provides it. On Windows,
`scripts/ensure_gram_native.jl` runs the supplied `build_gram.cmd` helper;
follow the Windows prerequisites for your GRAM distribution.

The wrapper checkout alone cannot build the library. Copy the NASA source and
runtime files first, as described above: the Unix helper needs
`Build/setup_cspice.sh` and the GRAM makefiles from that tree.

### Build the shared library

From the SpaceAGORA repository root, run:

```text
julia --project=. scripts/ensure_gram_native.jl
```

The script takes the GRAM root from an optional path argument, then from the
`GRAM_ROOT` environment variable, and otherwise uses
`data/GRAMSuite.jl/GRAM Suite 2.0`. After checking that the root and the
build helper exist, it skips the build if `Build/lib/libGRAM.<ext>` is already
there. It does not inspect that file's architecture or test whether it loads.
Otherwise it runs the build helper inside the GRAM tree:
`simulation/GRAM/build_gram.sh` on Unix or `simulation/GRAM/build_gram.cmd`
on Windows.

The Unix helper prepares CSPICE through `Build/setup_cspice.sh`, runs
`make shared` with the host's make command, and writes `gram.env` and
`.gram-build-manifest` next to itself. Both files hold local build paths; do
not commit them or copy them to another machine as installation instructions.
Use the helper rather than plain `make`: SpaceAGORA needs the shared library
that the `shared` target builds.

The result is:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib/libGRAM.so
```

On macOS the file ends in `.dylib`, on Windows in `.dll`. How long the build
takes depends on the host and the GRAM sources.

### Rebuild a copied or stale library

If the GRAM tree came from another machine with `Build/lib` already
populated, force a clean rebuild:

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

Use the same command when you suspect a stale local build after moving the
checkout. The Unix helper compares the host and root path stored at its last
build and rebuilds clean when either has changed, but `ensure_gram_native.jl`
does not run the helper at all when the library already exists and `--clean`
was not given. The message `Native GRAM library already present for this host`
is therefore not a compatibility check.

## Verify the asset layout

Run the asset report:

```text
julia --project=. src/interfaces/cli/main.jl assets check
```

On a machine set up for GRAM, the report shows both of these as `usable`:

- `gram_root`
- `spice_directory`

`missing-optional` for either means the submodule is not checked out, and
`pointer-only-optional` means its Git LFS files have not been fetched. The
report only checks that payload files are present. The wrapper alone already
makes both entries `usable` once its LFS files are fetched, so the report
cannot tell you whether the NASA files are in place, and it does not test
whether the native library loads. The native build above and the GRAM-backed
run below check those.

## First GRAM-backed run

Once the files are copied and the native library is built, try the basic
GRAM-backed example:

```text
julia --project=examples examples/AGORA_Basic_GRAMEarth.jl
```

It is the smallest SpaceAGORA run that exercises the `GRAMSuite.jl`
integration, the local GRAM runtime tree, the SPICE-backed Earth constructor
`Earth("", SPICE_PATH)` and `GRAMAtmosphereModel(planet_name="earth")`. For
larger mission examples, continue with `examples/AGORA_Vex.jl` (Venus) and
`examples/AGORA_Odyssey.jl` (Mars).

## GRAM-backed atmosphere

The keyword constructor of `GRAMAtmosphereModel` comes from the
`SpaceAGORAGRAMSuiteExt` package extension, which Julia loads once `GRAMSuite`
is imported alongside SpaceAGORA. `GRAMAtmosphereModel` is not exported from
the root `SpaceAGORA` module; reach it as
`SpaceAGORA.SimulationModel.GRAMAtmosphereModel`. The examples instead call
`setup_gram_example!()` from `examples/common.jl`, which loads `GRAMSuite` from
the submodule and defines `GRAMAtmosphereModel` in the calling module. In a
script that includes `examples/common.jl`, as the examples do:

```julia
setup_gram_example!()   # from examples/common.jl; loads GRAMSuite and the extension
density_model = GRAMAtmosphereModel(planet_name="earth")
```

### Opt-in density perturbations

Native GRAM supplies its mean density by default.
`SPACEAGORA_GRAM_DENSITY_PERTURBATION` selects a perturbed density instead:

- `step` advances a separate density random walk at each accepted solver step
  and holds its factor between steps.
- `pass` samples the walk along a drag-free prediction of each atmospheric
  pass and interpolates its factor.
- `naive_rhs` reads the walk inside the right-hand side. It is a diagnostic
  negative control, unsuitable for trajectory results, and the distributed RHS
  density service declines it.
- `off`, the default, leaves the density unperturbed.

`step` and `pass` apply the factor to spacecraft density on both the scalar
and the distributed density service paths, and only at or below the entry
interface. Temperature and wind are not perturbed by this setting.

Within one run, consecutive checkpoint segments keep the same perturbation
walk, pass counters and diagnostic history. Resuming from a saved checkpoint
is refused, because the checkpoint does not hold the native random state.
Multirate splitting is refused too, before its first subsolve, because its
overlapping intervals would revisit earlier walk times. Use another solver
mode or set the perturbation mode to `off`.

`SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_RESEED=1` derives each pass's
initial random stream from the seed and the pass index. The density history
still depends on what is sampled: the accepted steps in `step` mode and the
predicted trajectory in `pass` mode. In `pass` mode the last factor is held if
the real pass outlasts its prediction. These modeling choices and their
numerical accuracy need evidence for the intended scenario; the lifecycle and
interface tests do not validate a dispersed campaign.

### GRAM epoch alignment

A GRAM model carries its own epoch: the `initial_time` keyword of
`GRAMAtmosphereModel`, or the GRAMSuite default when it is omitted. The
simulation carries another in `SimulationConfiguration.initial_time`. Before
every run, `run_simulation` passes the density model through
[`with_density_model_epoch`](@ref) so the two agree. This happens before the
configuration is copied for state isolation, and the caller's configuration is
left unchanged.

For a `GRAMAtmosphereModel` built with keywords, an unchanged epoch returns
the same model object. The epochs are compared on their integer year, month,
day, hour and minute and on the seconds as `Float32`, with no time tolerance.
A changed epoch rebuilds a fresh native model from the recorded constructor
keywords with only `initial_time` replaced, so the planet, paths,
perturbation options and every other recorded keyword are kept. The rebuild
runs under the GRAM setup lock.

A `GRAMAtmosphereModel` wrapped around a raw GRAMSuite core has no recorded
constructor keywords and cannot be realigned; it raises an `ArgumentError`
when the epochs differ. Build it with the keyword constructor at the run's
`initial_time` instead.

A `GRAMAtmosphereModelSurrogate` built over a GRAM base keeps its object, file
and fallback setting for an unchanged epoch and rejects a changed one with an
`ArgumentError`: its table was generated for one epoch, and rebuilding the
native fallback would not re-epoch the table. Build and validate a surrogate
for the epoch you intend to run. Every other density model, including
surrogates over other base models, the grid models, `NRLMSISE00AtmosphereModel`
and the analytic models, passes through unchanged.

To avoid the rebuild, construct the model at the run's epoch. Continuing the
example above, with `t0` the `InitialTime` you also give the configuration as
`initial_time`:

```julia
t0 = InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)
density_model = GRAMAtmosphereModel(planet_name="earth", initial_time=t0)
```

The hook changes the construction epoch only. It does not convert between
time systems, validate atmospheric coordinate or datum conventions, certify
cached or surrogate data, or carry over a native random stream that was
advanced, or a handle that was edited by hand, before the run. The usual
native cache and environment settings apply to the rebuilt model as to any
other.

## Troubleshooting

### GRAM is reported missing or pointer-only

If `data/GRAMSuite.jl` exists but is empty, the submodule was never checked
out, and the asset report shows `gram_root` and `spice_directory` as
`missing-optional`. Run the `git submodule update` command above. Once the
wrapper is checked out, both entries are `pointer-only-optional` until its LFS
files are fetched and `usable` after that, whether or not the NASA files have
been copied. If the native build or a GRAM-backed example then fails, check
that the NASA GRAM Suite files were copied into
`data/GRAMSuite.jl/GRAM Suite 2.0`.

### `SPICE` is missing

The wrapper provides the `SPICE` folder, with its kernels in Git LFS. If the
asset report shows `spice_directory` as missing, check out the submodule; if
it shows it as pointer-only, run `git lfs pull` inside `data/GRAMSuite.jl`.
The NASA delivery has a `SPICE` folder of its own. Whatever you copy, the
folder must end up at:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
```

### The GRAM shared library is missing

Run:

```text
julia --project=. scripts/ensure_gram_native.jl
```

If the command reports that the library is present but the example still
cannot load it, one possible cause is a binary copied from another host. Try
the clean rebuild described above. If loading still fails, keep the complete
build or load error so the missing dependency or incompatible binary can be
identified.

### The build reports `gmake not found` on macOS

Install GNU Make, then rerun the build command:

```text
brew install make
```

### The native build files are missing

If the build cannot find `Build/setup_cspice.sh` or the GRAM makefiles, check
that the NASA GRAM delivery was copied to the expected location. The wrapper
and its SPICE kernels alone are not enough to compile `libGRAM`.

### The GRAM package extension fails to load

Treat an `Error during loading of extension SpaceAGORAGRAMSuiteExt` message as
an incomplete setup, even if the script carries on. Keep the complete error
and the SpaceAGORA and GRAMSuite revision identifiers, and check that the two
revisions are supported together before you use the run's results: the
supported pairing is the submodule commit that SpaceAGORA records. A failed
extension can leave atmosphere or ephemeris hooks unavailable. Updating
GRAMSuite to an arbitrary branch tip does not establish compatibility.

### Git LFS reports `no space left on device`

The submodule holds large LFS files. Free disk space first, then retry the LFS
checkout:

```text
cd data/GRAMSuite.jl
git lfs pull
git checkout .
```

If the failed checkout left temporary LFS files behind, you can remove them
once you are sure no other Git LFS operation is running. These paths assume an
ordinary clone, where `.git` at the repository root is a directory:

```text
rm -rf ../../.git/modules/GRAMSuite.jl/lfs/incomplete/*
rm -rf ../../.git/modules/GRAMSuite.jl/lfs/tmp/*
```
