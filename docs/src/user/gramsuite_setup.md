# GRAMSuite Setup

Use this page when you have been granted access to NASA GRAM and need to wire
the licensed files into a local SpaceAGORA checkout.

This page is for users running higher-fidelity GRAM/SPICE-backed workflows. It
is not required for the baseline no-GRAM onboarding path.

What to read next:

- [Installation & Environment](installation_environment.md)
- [Assets & Modes](../assets.md)
- [Verification Study](verification_study.md)

## What this document covers

This guide explains how to:

1. pull the vendored `GRAMSuite.jl` submodule used by SpaceAGORA
2. request and receive the official NASA GRAM Suite distribution
3. copy the GRAM installation folders into the repository-local expected path
4. verify that SpaceAGORA can see the resulting GRAM and SPICE directories

It does not distribute GRAM itself. Access to the official GRAM Suite must be
requested separately from NASA.

## Request GRAM access

Request access through the NASA Software Catalog entry for GRAM:

- <https://software.nasa.gov/software/MFS-33888-1>

After approval, obtain the official GRAM Suite package through the process
specified by NASA. Keep the original distribution intact until you finish the
copy/verification steps below.

## Pull the `GRAMSuite.jl` submodule

Start from the SpaceAGORA repository root and make sure the submodule is
present locally. The full checkout uses Git LFS and downloads several GB of
binary GRAM data. Before running it, make sure the machine has enough free disk
space for both the final files and Git LFS temporary downloads; 15-20 GB free is
a practical minimum.

```text
git submodule update --init --recursive --remote
```

This should populate the vendored Julia wrapper path:

```text
data/GRAMSuite.jl
```

The wrapper package is not the same thing as the licensed GRAM binaries and
data. The submodule gives SpaceAGORA the Julia integration layer and the
expected folder scaffold; the official NASA distribution provides the GRAM
content that must be placed into that scaffold.

If you only need the wrapper source or want to inspect the scaffold without
downloading LFS objects immediately, skip LFS smudging during the submodule
checkout:

```text
GIT_LFS_SKIP_SMUDGE=1 git submodule update --init --recursive --remote
```

Later, after freeing enough disk space, fetch the LFS-backed GRAM files from
inside the submodule:

```text
cd data/GRAMSuite.jl
git lfs pull
```

## Expected target location

SpaceAGORA expects the official GRAM Suite tree under:

```text
data/GRAMSuite.jl/GRAM Suite 2.0
```

The most important subpaths for SpaceAGORA are:

- `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE`
- `data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib`
- the rest of the official GRAM Suite runtime tree required by `GRAMSuite.jl`

## Copy the official GRAM Suite folders

Once NASA has provided access and you have the official GRAM Suite files
available locally, copy the contents of the official distribution into the
vendored target location inside this repository.

At a high level, the process is:

1. unpack or open the official GRAM Suite delivery from NASA
2. locate the top-level `GRAM Suite 2.0` directory from that delivery
3. copy the contents of that directory into `data/GRAMSuite.jl/GRAM Suite 2.0`
4. confirm that the final path is exactly
   `data/GRAMSuite.jl/GRAM Suite 2.0`

After the copy, these should exist:

```text
data/GRAMSuite.jl/GRAM Suite 2.0
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
data/GRAMSuite.jl/GRAM Suite 2.0/simulation
```

If the official delivery already contains a `GRAM Suite 2.0` directory, copy
that directory as-is. Avoid renaming it, because SpaceAGORA and the wrapper
scripts expect that exact folder name.

## Build the native GRAM library

GRAM is a C/C++ codebase that SpaceAGORA calls through FFI. The shared library
is compiled per machine and is never shipped prebuilt, so this step is required
on every host.

### Build prerequisites

| Requirement | Notes |
|---|---|
| C/C++ toolchain | GCC on Linux; Clang via Xcode command-line tools on macOS; MinGW-w64 (MSYS2) on Windows |
| GNU Make | Linux: `make`. macOS: `gmake`, from `brew install make` — the build files use GNU Make features that `/usr/bin/make` does not support. Windows: `mingw32-make`. |
| CSPICE | Bundled and selected automatically on Linux x86_64 and Windows MinGW. On macOS, run `brew install cspice` first. On Linux arm64, supply your own archive (see below). |

A Fortran compiler is **not** required, despite `gfortran` appearing in the
GRAM build configuration. `Build/makefile.defs` ends with `undefine FC`, which
disables the Fortran example targets; the library build never invokes it.

### Build it

From the repo root:

```text
julia --project=. scripts/ensure_gram_native.jl
```

This is the only command you need. It returns immediately when the library for
your host already exists, so it is safe to run before every GRAM-backed job.
Otherwise it delegates to the vendored helper
`GRAM Suite 2.0/simulation/GRAM/build_gram.sh`, which:

1. runs `Build/setup_cspice.sh` to put the right CSPICE archive in place —
   you do not run this yourself
2. runs `make shared` with one job per core, using the correct make binary for
   the host
3. writes the library to
   `data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib/libGRAM.so` (`.dylib` on macOS,
   `.dll` on Windows)
4. writes `simulation/GRAM/gram.env` and `simulation/GRAM/.gram-build-manifest`,
   recording the host and root path the artifacts belong to

A full build takes roughly 20 seconds on 24 cores, or about 3 minutes of total
CPU time.

`shared` is the only make target that produces `libGRAM`. Plain `make` and the
per-planet targets such as `make Mars` build static libraries and example
executables but leave `Build/lib` without the shared library SpaceAGORA loads.
There is no faster single-planet path to a working `libGRAM`.

### When you need `--clean`

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

Use it when **the GRAM tree was copied or `rsync`ed in from another machine
with `Build/lib` already populated**. `scripts/remote/spaceagora-remote`
mirrors the whole working tree, so remote hosts hit this routinely. The
ordinary command's first test is whether `Build/lib/libGRAM.<ext>` exists, and
a library built on a different machine satisfies that test — the build is
skipped and you get a load error, or worse, a binary built for the wrong
architecture.

You do **not** need `--clean` for a checkout that moved or was renamed on the
same machine. The build manifest records the host tag and the root path, and
the helper forces a clean rebuild by itself when either has changed, printing
what it detected.

### CSPICE on Linux arm64

No archive is bundled for arm64/aarch64. Either install CSPICE so that one of
`/usr/lib/libcspice.a`, `/usr/local/lib/libcspice.a`, `/usr/lib64/libcspice.a`,
or `/usr/lib/aarch64-linux-gnu/libcspice.a` exists, or point the build at an
explicit archive:

```text
cd "data/GRAMSuite.jl/GRAM Suite 2.0/Build"
make shared -j SPICE_LIB=/absolute/path/to/cspice.a
```

## Verify the asset layout

Use the built-in asset report:

```text
julia --project=. src/cli/main.jl assets check
```

For a GRAM-ready machine, the report should show these as available:

- `gram_root`
- `spice_directory`

If either is still missing, re-check the final directory names and nesting
under `data/GRAMSuite.jl`.

## First GRAM-backed run

Once the folder copy and native-library setup are complete, try the basic
GRAM-backed example:

```text
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

That example is the smallest SpaceAGORA path that exercises:

- the vendored `GRAMSuite.jl` integration
- the local GRAM runtime tree
- the SPICE-backed Earth constructor
- `GRAMAtmosphereModel(planet_name="earth")`

If you want a larger mission example after that, move on to:

- `examples/AGORA_Vex.jl`
- `examples/AGORA_Odyssey.jl`

## Troubleshooting

### The submodule exists, but GRAM is still reported missing

This usually means only the `GRAMSuite.jl` wrapper was pulled, but the official
NASA GRAM Suite folders were not copied into
`data/GRAMSuite.jl/GRAM Suite 2.0`.

### `SPICE` is missing

Make sure the official GRAM delivery's `SPICE` folder ended up at:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/SPICE
```

### The GRAM shared library is missing

Run:

```text
julia --project=. scripts/ensure_gram_native.jl
```

If it prints `Native GRAM library already present for this host` but loading
still fails, the library on disk was built somewhere else. Force a rebuild:

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

### The build fails with `Could not auto-find a Linux CSPICE archive`

You are on an architecture with no bundled CSPICE. See
[CSPICE on Linux arm64](#cspice-on-linux-arm64).

### The build fails with `gmake not found` on macOS

```text
brew install make
```

The GRAM makefiles require GNU Make; the `make` shipped with macOS is too old.

### `make` succeeded but there is still no `libGRAM`

You most likely ran `make` or a per-planet target rather than `make shared`.
Only the `shared` target links the FFI library. Prefer
`scripts/ensure_gram_native.jl`, which always builds the right target and also
records the build manifest that the staleness detection depends on.

### Git LFS reports `no space left on device`

The GRAM submodule contains large LFS-backed binaries. Free disk space first,
then retry the LFS checkout:

```text
cd data/GRAMSuite.jl
git lfs pull
git checkout .
```

If the failed checkout left temporary LFS files behind, they can be removed
after confirming no other Git LFS operation is running:

```text
rm -rf ../../.git/modules/GRAMSuite.jl/lfs/incomplete/*
rm -rf ../../.git/modules/GRAMSuite.jl/lfs/tmp/*
```
