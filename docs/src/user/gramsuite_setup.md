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

The submodule tracks the public `GRAMSuite.jl` mirror, which needs no special
access. It carries the Julia integration layer, `simulation/`, `SPICE/` and the
folder scaffold — but not the licensed GRAM content. The official NASA
distribution provides that, and it must be placed into the scaffold; see
[Copy the official GRAM Suite folders](#copy-the-official-gram-suite-folders).

### Lab members: using `dev-GRAMSuite.jl` instead

`dev-GRAMSuite.jl` carries the complete GRAM Suite tree — `Build/`, `common/`,
every planet directory — so pointing at it removes the copy step entirely. It is
private, and requires access to the GRAM code.

Override the URL **for your machine only**. Do not commit a change to
`.gitmodules`: the public URL has to keep resolving for everyone without dev
access, who would otherwise be unable to initialize the submodule at all.

```text
DEV=https://github.com/Space-FALCON-Lab/dev-GRAMSuite.jl.git
git config submodule.GRAMSuite.jl.url "$DEV"
git -C data/GRAMSuite.jl remote set-url origin "$DEV"
git -C data/GRAMSuite.jl fetch origin main
git -C data/GRAMSuite.jl checkout origin/main
```

Both `config` lines are needed: the first is what `git submodule` commands
consult, the second is what a `git -C data/GRAMSuite.jl fetch` uses.

Three consequences, all normal for this setup:

**Do not run `git submodule sync` afterwards.** That command copies the URL
*from* `.gitmodules` *into* `.git/config`, so it silently puts you back on the
public mirror. Its legitimate use is the opposite case, below.

**`git submodule update` reverts the checkout.** The recorded gitlink names a
public commit, so a plain `git submodule update` — including one inside a setup
script or CI step — checks the public content back out. Re-run the
`fetch`/`checkout` pair above if that happens.

**`git status` shows `data/GRAMSuite.jl` as modified, permanently.** The
submodule sits on a commit that is not the recorded one. Do not commit that
pointer: the two repositories share no objects, so a dev commit is never a valid
gitlink for consumers resolving against the public URL, and committing one
breaks every clone without dev access.

### When `git submodule sync` *is* needed

If the URL in `.gitmodules` changes upstream and you already have an initialized
checkout, pulling is not enough on its own. `git submodule init` copied the old
URL into `.git/config` at setup time, and that copy is what fetches use.
Propagate it with:

```text
git submodule sync
git submodule update --init --recursive
```

A fresh clone reads `.gitmodules` directly and needs neither command.

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

## Verified end-to-end walkthrough

The sequence below was run start to finish on a clean machine
(Ubuntu 24.04, x86_64, 24 cores, Julia 1.12.1, GCC 13.3.0) from an empty
directory. Timings are from that run and are indicative, not guarantees.

```text
# 1. clone and instantiate                                            ~40 s
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"

# 2. confirm the baseline path works before involving GRAM at all     ~25 s
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --smoke

# 3. fetch the GRAMSuite wrapper (Git LFS; ~3 GB, several minutes)    ~3 min
git submodule update --init --recursive --remote

# 4. copy the licensed NASA GRAM Suite delivery into the scaffold     ~40 s
#    (see "Copy the official GRAM Suite folders" above)
#    the tree is ~17 GB once complete

# 5. build the native library                                         ~21 s
julia --project=. scripts/ensure_gram_native.jl

# 6. confirm the assets are visible
julia --project=. src/cli/main.jl assets check

# 7. first GRAM-backed run                                            ~34 s
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

Two points that are easy to misread:

Step 3 does **not** give you a buildable GRAM tree, and step 4 is not optional
on the default path. The public mirror the submodule tracks strips the
proprietary folders — no `Build/`, no `common/`, no planet directories — so
running step 5 before step 4 fails with nothing to build. Lab members who point
the submodule at `dev-GRAMSuite.jl` get the complete tree and can skip step 4.

That is a property of the public mirror, not of GRAMSuite itself. Lab members
can point the submodule at `dev-GRAMSuite.jl`, which carries the complete tree,
and skip step 4 — see
[Lab members: using `dev-GRAMSuite.jl` instead](#lab-members-using-dev-gramsuitejl-instead).

Step 5 is the only step that compiles anything, and it is fast. If it runs for
minutes, or if it reports a path that is not inside your checkout, stop and read
the troubleshooting entries below rather than waiting.

Expected output from step 6 on a GRAM-ready machine:

```text
- gram_root: available
- spice_directory: available
```

Step 7 writes `output/simulation_results.csv`. A successful run prints
`COMPUTATIONAL TIME` and no `Error during loading of extension` line — see
[The GRAM extension fails to load](#the-gram-extension-fails-to-load) if it does.

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

### The GRAM extension fails to load

**Symptom:** the run completes, but begins with

```text
Error: Error during loading of extension SpaceAGORAGRAMSuiteExt of SpaceAGORA
InitError: UndefVarError: `_GRAM_EPHEMERIS_STATE_FN` not defined in `GRAMSuite`
```

**Cause:** the `GRAMSuite` checkout is older than the extension in this
repository. The extension's `__init__` installs several hooks and aborts on the
missing one, so the hooks after it — including the process-wide CSPICE lock that
GRAM model construction must serialise against — are never installed. Because
this is reported as a warning rather than an error, the simulation continues
without them. Single-threaded results are unaffected; what is lost is thread
safety, so a threaded or multi-satellite GRAM run can corrupt SPICE state.

**Resolution:** update the submodule to a `GRAMSuite` revision that defines the
symbol named in the error:

```text
cd data/GRAMSuite.jl
git fetch origin
git checkout origin/main
```

Never ignore this warning on a threaded run.

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
