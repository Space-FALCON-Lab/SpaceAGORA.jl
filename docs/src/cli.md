# CLI

The command-line interface runs examples, verification studies and benchmark
launchers, opens the viewer on saved results, and inspects assets, without
writing a Julia script. Its entry point is a Julia file, so the same command
works on Windows, Linux and macOS:

```text
julia --project=. src/interfaces/cli/main.jl assets check
```

Two wrappers start it for you: `./bin/spaceagora` on Linux and macOS, and
`bin\spaceagora.bat` on Windows. Each starts Julia with `--project` set to
`JULIA_PROJECT`, or to the repository root when `JULIA_PROJECT` is unset.
`help` (or `--help`, `-h`, or no arguments at all) prints the usage summary.
From Julia, `SpaceAGORA.run_cli(args)` runs the same commands and returns the
exit status.

See also the [Quickstart](user/quickstart.md), the [Recipes](user/recipes.md)
and [Assets & Modes](assets.md).

## What each command needs

`run`, `telemetry` and `benchmark` start a child Julia process. Unless a command
says otherwise below, the child runs in the examples environment
(`examples/Project.toml`), which provides this checkout's `SpaceAGORA` and the
packages the examples and studies use beyond it, such as Plots. Prepare it once
with `julia --startup-file=no examples/setup.jl`
([Installation & Environment](user/installation_environment.md#Examples-environment));
until you do, these commands stop with a message that names that command.

The rest depends on the script the command starts, not on the CLI:

| Command | Needs on top of the instantiated repository and the examples environment |
|---|---|
| `assets check`, `assets manifest`, `assets setup-open` | nothing |
| `run --example=<no-GRAM example>` (`AGORA_Basic_Quickstart.jl`, `AGORA_Earth_NoGRAM.jl`, `AGORA_Earth_MonteCarlo.jl`, `Solar_Panel_Cloth_Deployment_Demo.jl`) | nothing |
| `run --example=<GRAM-backed or SPICE-backed example>` (`AGORA_Earth_Aerobraking.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl`, `Earth_Thruster_Test.jl`, `AGORA_Keplerian.jl`) | the `data/GRAMSuite.jl` submodule ([GRAMSuite Setup](user/gramsuite_setup.md)); the GRAM-backed examples also need the native GRAM library built |
| `run --example=` an EDG example (`odyssey_surrogate.jl`, `AGORA_Odyssey_Energy_Depletion_Control_Test.jl`) | a prepared project that declares `SpaceAGORA` and `SpaceAGORAEDG`; see [Run an example](#Run-an-example) |
| HYPR and RPO examples | the prepared `examples/rpo_planner_env` environment, and SPICE kernels for the RPO cases. Run them directly with the commands in the [Examples Catalog](user/examples_catalog.md#RPO-and-robotics); the CLI does not select that environment. |
| `telemetry ...` | the `data/GRAMSuite.jl` submodule: the study loads the vendored `GRAMSuite` package before it reads any scenario, even with `--scenarios=odyssey`. The truth files it grades against are in the repository. |
| `benchmark ...` | the `data/GRAMSuite.jl` submodule, and the native GRAM library for the GRAM-backed cases |

Without the submodule, `telemetry` and `benchmark` stop with "Package
GRAMSuite not found in current path". That is the missing prerequisite, not a
CLI fault.

## Commands

### Run an example

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

| Option | Effect |
|---|---|
| `--example=<file>` | an example file name under `examples/` (the `.jl` may be left off), or a path to a script |
| `--output-dir=<dir>` | sets `SPACEAGORA_CLI_OUTPUT_DIR` for the child process |
| `--smoke` | sets `SPACEAGORA_EXAMPLE_SMOKE=1` and `SPACEAGORA_EXAMPLE_SMOKE_RESULTS=1`, the shortened mode that long examples support |
| `--visualize` | sets `SPACEAGORA_VISUALIZATION=1`, loads the viewer in the selected project and exports a page after the run |
| `--print-only` | prints the selected project, the script, the environment variables and the full command, then exits without running |

Without `--visualize`, the child runs in the examples environment, even if the
CLI itself was started with another project. The two EDG examples are the
exception: they run in the project named by `SPACEAGORA_EDG_ENV`, or else in the
project the CLI was started with, and that project must declare both
`SpaceAGORA` and `SpaceAGORAEDG` (with `--visualize`, the viewer project
described below must declare `SpaceAGORAEDG` as well). The prepared
`examples/odyssey_surrogate_env` from the
[Odyssey tutorial](tutorials/odyssey_surrogate.md) qualifies:

```text
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl run --example=odyssey_surrogate.jl
```

### View saved results or export a page after a run

The viewer lives in a separate project. Prepare it once, then point the CLI at
it:

```sh
julia --startup-file=no scripts/setup_visualizer.jl ../spaceagora-viewer
julia --project=../spaceagora-viewer src/interfaces/cli/main.jl visualize --run=output/cli_run
julia --project=../spaceagora-viewer src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_viewer --visualize
```

`visualize` reads a saved result bundle and its scene without running the
simulation again. `run --visualize` asks the example to record its scene and
export the HTML page. The entry point loads `SpaceAGORAVisualizer` for these two
commands only; when you call `run_cli` from Julia to view saved results, load
`SpaceAGORA` and `SpaceAGORAVisualizer` yourself first.

| `visualize` option | Effect |
|---|---|
| `--run=<path>` | required: a results directory holding `simulation_results.feather`, the bundle prefix itself, or with `--ensemble` a campaign directory of sample subdirectories |
| `--out=<html>` | where to write the page |
| `--max-frames=<n>` | the largest number of frames the page holds |
| `--frame=inertial` or `--frame=planet_fixed` | the frame the scene is drawn in |
| `--texture=best`, `4k` or `8k` | the planet texture resolution; `--no-textures` turns textures off |
| `--trail-orbits=<n>` | trail length, in orbits |
| `--title=<text>` | the page title |
| `--model=<id>=<file>` | an STL, OBJ or GLB model for spacecraft `<id>`; repeat for more spacecraft. Not with `--ensemble`. |
| `--model-scale=<m per unit>`, `--model-center=0` or `1` | the model's scale, and whether to center it on the spacecraft |
| `--ensemble` | build the ensemble page for a campaign directory |

The command prints `viewer=<path>` for the page it wrote.
[Visualization](user/visualization.md) describes the saved payload and the
viewer options in detail.

A visualized `run` child uses the project the CLI was started with, unless
`SPACEAGORA_VISUALIZER_ENV` names another prepared project. Either way the
project must declare both `SpaceAGORA` and `SpaceAGORAVisualizer`; the CLI
checks this before it launches anything, `--print-only` included. The child
loads the viewer before the example starts, so an example that activates the
core project keeps viewer support, and the script's `PROGRAM_FILE`, arguments
and main guard behave as usual. If the project declares the packages but they
are not installed, the preload fails before the example writes any output. The
CLI never installs packages. With the wrappers, which start under
`JULIA_PROJECT` or the repository root, set `SPACEAGORA_VISUALIZER_ENV` (or
`JULIA_PROJECT`) to the prepared viewer project.

### Run telemetry verification

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
julia --project=. src/interfaces/cli/main.jl telemetry --profile=full --output-dir=output/telemetry_full
```

The profile is `quick` (the default), `full`, or `smoke`, which is another name
for `quick`. Give it as a bare word or as `--profile=<name>`.

| Option | Effect |
|---|---|
| `--output-dir=<dir>` | where the CSV files go (default `output/telemetry`) |
| `--enforce=0` or `1` | with `1`, the command fails when a verification threshold fails |
| `--plots=0` or `1` | turns plot generation on or off; off by default |
| `--scenarios=a,b` | runs only the named manifest scenarios (default: all); an unknown name is an error |
| `--print-only` | prints the resolved command without running it |

Plots are off by default so that the command still works where the plotting
scripts are not available. The study writes `telemetry_orbit_accuracy_summary.csv`
and `telemetry_orbit_accuracy_errors.csv` to the output directory.

### Run benchmark launchers

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_cli
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
```

| Mode | Launcher | Output directory variable (default) |
|---|---|---|
| `runtime-analysis` | `benchmarks/studies/performance_runtime_analysis.jl` | `SPACEAGORA_PERF_OUTDIR` (`output/performance`) |
| `smart-parallel-ladder` | `benchmarks/studies/performance_smart_parallel_ladder.jl` | `SPACEAGORA_SMART_LADDER_OUTDIR` (`output/performance/smart_parallel_ladder`) |

`--output-dir=<dir>` sets that variable and `--print-only` prints the command
without running it. Every other argument after the mode, such as the profile
`quick`, `full` or `smoke`, goes to the launcher unchanged:

```text
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder full --passes=3 --output-dir=output/smart_ladder_full
```

### Inspect assets

```text
julia --project=. src/interfaces/cli/main.jl assets check
julia --project=. src/interfaces/cli/main.jl assets manifest
julia --project=. src/interfaces/cli/main.jl assets setup-open
```

`assets check` reports whether each asset root and catalog in this checkout is
usable, missing, empty, an undownloaded Git LFS pointer or empty of payload,
and exits with status 1 when a required entry is not usable. It confirms that
payload bytes exist, not that a grid is valid. `assets manifest` prints
`data/assets_manifest.toml`, and `assets setup-open` reports the open-asset
setup. [Assets & Modes](assets.md) explains the states.

## Print-only mode

`run`, `telemetry` and `benchmark` accept `--print-only`, which prints the
project, the script, the environment variables the CLI would set and the exact
command, on lines starting `project=`, `script=`, `env:` and `cmd=`, and runs
nothing:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --smoke --print-only
julia --project=. src/interfaces/cli/main.jl telemetry smoke --output-dir=output/telemetry_smoke --print-only
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke --print-only
```

### Named surrogate data

After [setting up the Odyssey environment](tutorials/odyssey_surrogate.md):

```sh
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets list
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets fetch --preset=odyssey_p20_frozen_v1 --version=1.0.0
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets check --preset=odyssey_p20_frozen_v1 --version=1.0.0
```

`list` reads the catalog without downloading anything. `fetch` downloads only
the selected version and verifies its bytes; `fetch --offline` verifies an
installed copy without network access. `check` works offline on installed data
and validates the payload's schema, coordinates or domain limits and generation
settings. `--file=<path>` selects an explicit copy, which must match the
preset byte for byte; an invalid path is an error and never falls back to a
download. Options can be written `--preset=<id>` or `--preset <id>`, and
likewise for `--version` and `--file`. An exact version is always required.

The same commands take `--preset=mars_global_upper_p20_frozen_v1 --version=1.0.0`
for the global upper atmosphere and
`--preset=mars_global_near_surface_p20_frozen_v1 --version=1.2.0` (or
`--version=1.1.0` or `--version=1.0.0`) for the near-surface atmosphere.
`--preset=mars_global_descent_p20_frozen_v1 --version=1.0.0` selects the
composite descent preset, which has no payload of its own: `fetch` installs and
verifies both components, `check` builds the composite from them, and neither
takes `--file`.

A preset `check` loads the public `GRAMSuite` package (the Julia wrapper, not
native GRAM) and refuses to run in a project that cannot load it; neither
command builds or initializes native GRAM. The Odyssey example retrieves its
preset itself when it starts, so these commands are only for inspecting or
prefetching. The plain `assets check` above covers the repository's own asset
roots, not preset contents. The `runtime_wrapper_revision` that `check` prints
is the GRAMSuite revision the preset was generated with; the revision your
environment runs is pinned in its `Project.toml`.
