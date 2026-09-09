# Quickstart

Use this page if you want the fastest path from clone to a successful local
run.

This page is for first-time users who do not need GRAM or SPICE on day one.

This default quickstart baseline is the Earth vacuum path:
`NoAtmosphereModel()` plus `SimpleEphemeridesModel()`.

Shortest successful path:

```text
GIT_LFS_SKIP_SMUDGE=1 git clone --filter=blob:none https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Assets & Modes](../assets.md)
- [First Simulation](first_simulation.md)
- [Simulation Outputs](outputs.md)
- [Recipes](recipes.md)

## What to expect

Measured on a fresh clone (September 2026, packages already in the local
depot): the clone takes a few minutes, `Pkg.instantiate()` about a minute, and
the first run about a minute, most of it compilation. The run prints the
initial state, the number of saved samples and the computational time, and
writes:

```text
output/
  simulation_results.csv
  simulation_results.feather
  simulation_results.manifest.toml
  plots/
    quickstart_altitude_speed.png
    quickstart_inertial_trajectory.png
    quickstart_3d_orbit.png
    quickstart_ground_track.png
```

`output/` is the repository root's `output/` directory (gitignored). Every
example writes the same three file names there, so a second run overwrites
the first; see [Simulation Outputs](outputs.md) for keeping runs apart.

`GIT_LFS_SKIP_SMUDGE=1` on the clone line keeps the seven GRAM surrogate grids
under `data/GRAM_surrogate/` as small Git LFS pointers. They are about 2.5 GB,
nothing on this path reads them, and `git lfs pull --include
"data/GRAM_surrogate/*"` fetches them later if a benchmark study needs them.

## What this run gives you

The no-GRAM example is the shortest repository-owned path to a working local
simulation. It does not require:

- a local GRAM installation
- SPICE kernels
- surrogate atmosphere grids
- the `data/GRAMSuite.jl` submodule (it stays empty in a fresh clone)

## If you prefer the CLI

The equivalent command-line path is:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

`--output-dir` is what keeps this run's files apart from the script run above:
they land under `output/cli_run/` with the same names. Add `--smoke` for a
two-minute mission that only checks the path end to end; it also writes to
`output/` (or the directory you give) and overwrites what is there.

## When to stop using this page

Move on as soon as one of these becomes true:

- you need higher-fidelity assets or licensed data
- you want to run a fuller scenario than the baseline no-GRAM example
- you want to inspect the files written under `output/`
- you want the verification study rather than a first smoke run
