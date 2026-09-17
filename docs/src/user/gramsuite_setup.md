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

## Build or verify the native GRAM library

SpaceAGORA loads GRAM through a native shared library, `libGRAM`. A library
copied from another machine may not work on this host, even when its filename
has the expected extension.

### Build prerequisites

Use the compiler and CSPICE archive appropriate for your licensed GRAM
distribution and host architecture. The Unix build helper uses GNU Make:
`make` on Linux and `gmake` on macOS. On macOS, install it with
`brew install make` if needed. Windows uses the supplied `build_gram.cmd`
helper; follow the Windows prerequisites for your GRAM distribution.

The wrapper checkout alone is not a complete native build environment. Copy
the licensed source and runtime files first, as described above. The Unix
helper expects `Build/setup_cspice.sh` and the GRAM makefiles in that tree.

### Build the shared library

From the SpaceAGORA repository root, run:

```text
julia --project=. scripts/ensure_gram_native.jl
```

After validating the GRAM root and build-helper path, this command skips the
build if `Build/lib/libGRAM.<ext>` already exists. It does not inspect the
binary's architecture or test whether it can load.
Otherwise it invokes the build helper inside the selected GRAM tree:
`simulation/GRAM/build_gram.sh` on Unix or `simulation/GRAM/build_gram.cmd`
on Windows.

The Unix helper prepares CSPICE through `Build/setup_cspice.sh`, runs
`make shared` using the host's make command, and writes `gram.env` and
`.gram-build-manifest` beside the helper. Those files contain local build paths
and should not be committed or copied as installation instructions for another
machine. Prefer the helper to running plain `make`: SpaceAGORA needs the shared
library that the `shared` target produces.

The expected result is:

```text
data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib/libGRAM.so
```

The filename ends in `.dylib` on macOS and `.dll` on Windows. There is no fixed
build time; it depends on the host and the supplied GRAM sources.

### Rebuild a copied or stale library

If the GRAM tree was copied from another machine with `Build/lib` populated,
force a clean rebuild:

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

The same command is useful when a stale local build is suspected after moving
the checkout. The Unix helper checks its stored host and root path when it
runs, but the Julia wrapper skips that helper whenever the expected library
already exists and `--clean` was not requested. Do not treat the
`Native GRAM library already present for this host` message as a compatibility
check.

## Verify the asset layout

Use the built-in asset report:

```text
julia --project=. src/cli/main.jl assets check
```

For a GRAM-ready machine, the report should show these as available:

- `gram_root`
- `spice_directory`

If either is still missing, re-check the final directory names and nesting
under `data/GRAMSuite.jl`. This report checks asset presence, not whether the
native library can load; verify that with the GRAM-backed run below.

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

If the command reports that the library is present but the example still
cannot load it, a binary copied from another host is one possible cause. Try
the clean rebuild described above. If it still fails, keep the complete build
or load error so the missing dependency or incompatible binary can be checked.

### The build reports `gmake not found` on macOS

Install GNU Make, then rerun the build command:

```text
brew install make
```

### The native build files are missing

If the build cannot find `Build/setup_cspice.sh` or the GRAM makefiles, check
that the licensed GRAM delivery was copied to the expected location. Having
the Julia wrapper and SPICE kernels alone is not enough to compile `libGRAM`.

### The GRAM package extension fails to load

Treat an `Error during loading of extension SpaceAGORAGRAMSuiteExt` message as
an incomplete setup even if the script continues. Keep the complete error and
the SpaceAGORA and GRAMSuite revision identifiers. Check that the two revisions
are supported together before using the run's results; a failed extension can
leave atmosphere or ephemeris hooks unavailable. Updating GRAMSuite to an
arbitrary branch tip does not establish compatibility.

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
