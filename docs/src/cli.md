# CLI

Use this page when you want the stable command-line surface instead of calling
Julia scripts directly.

This page is for users and operators who prefer packaged commands for examples,
verification studies, benchmarks, and asset inspection.

Shortest successful command:

```text
julia --project=. src/interfaces/cli/main.jl assets check
```

What to read next:

- [Quickstart](user/quickstart.md)
- [Recipes](user/recipes.md)
- [Assets and Modes](assets.md)

SpaceAGORA exposes a package-owned CLI surface through:

```text
julia --project=. src/interfaces/cli/main.jl
```

This Julia entrypoint is the cross-platform command path for Windows, Linux,
and macOS.

Convenience wrappers:

- Linux/macOS: `./bin/spaceagora`
- Windows: `bin\spaceagora.bat`

## What each command needs

The CLI normally starts child processes (`run`, `telemetry`, `benchmark`) in
the examples environment (`examples/Project.toml`), which provides this
checkout's `SpaceAGORA` and the packages the examples and studies use beyond it,
such as Plots. Prepare it once with `julia --startup-file=no examples/setup.jl`
([Installation & Environment](user/installation_environment.md#Examples-environment));
until then those commands stop and print that command. A `run --visualize`
request instead uses the prepared visualizer project described below, and an
example started there continues in the examples environment through
`examples/common.jl`. `--print-only` prints the selected project on its
`project=` line and the full command on its `cmd=` line without running
anything.

What a command needs beyond `Pkg.instantiate()` and the examples setup depends
on the script it starts, not on the CLI:

| Command | Needs on top of the instantiated repository and the examples environment |
|---|---|
| `assets check`, `assets manifest`, `assets setup-open` | nothing |
| `run --example=<no-GRAM example>` (`AGORA_Basic_Quickstart.jl`, `AGORA_Earth_NoGRAM.jl`, `AGORA_Earth_MonteCarlo.jl`, `Solar_Panel_Cloth_Deployment_Demo.jl`) | nothing |
| `run --example=<GRAM-backed or SPICE-backed example>` (`AGORA_Earth_Aerobraking.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl`, `Earth_Thruster_Test.jl`, `AGORA_Keplerian.jl`) | the `data/GRAMSuite.jl` submodule ([GRAMSuite Setup](user/gramsuite_setup.md)); the GRAM-backed ones also need the native GRAM library built |
| HYPR/RPO examples (direct script commands) | the prepared `examples/rpo_planner_env` companion environment; RPO cases also need SPICE kernels. Use the setup and direct commands in the [Examples Catalog](user/examples_catalog.md#RPO-and-robotics). |
| `telemetry ...` | the `data/GRAMSuite.jl` submodule: the study loads the vendored `GRAMSuite` package before it reads any scenario, even for `--scenarios=odyssey`; the truth files it grades are in the repository |
| `benchmark ...` | the `data/GRAMSuite.jl` submodule, and for the GRAM-backed cases the native GRAM library |

The CLI does not automatically select the HYPR companion environment. An
ordinary `run` launches the examples environment, even if the CLI itself was
started with `--project=examples/rpo_planner_env`. Run HYPR examples directly
with that companion project as shown in the catalog.

Without the submodule, `telemetry` and `benchmark` stop with "Package
GRAMSuite not found in current path"; that is the prerequisite, not the CLI.

## Commands

### Run an example

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

Supported options:

| Option | Effect |
|---|---|
| `--example=<file>` | Example filename under `examples/`, or a direct path |
| `--output-dir=<dir>` | Sets `SPACEAGORA_CLI_OUTPUT_DIR` for the child process |
| `--smoke` | Sets the example smoke-mode environment variables |
| `--visualize` | Loads the optional renderer in the selected project and exports a page after the run |
| `--print-only` | Prints the resolved launcher without running it |

Use `--smoke` for long examples that support
`SPACEAGORA_EXAMPLE_SMOKE=1`.

Without `--visualize`, the child runs in the examples environment.
`--print-only` shows the selected project and the exact command.

### View saved results or export a page after a run

Prepare the optional renderer once in a separate project:

```sh
julia --startup-file=no scripts/setup_visualizer.jl ../spaceagora-viewer
julia --project=../spaceagora-viewer src/interfaces/cli/main.jl visualize --run=output/cli_run
julia --project=../spaceagora-viewer src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_viewer --visualize
```

`visualize` reads a saved scene and result bundle without repeating the
simulation. `run --visualize` requests scene recording and HTML export from the
example. The direct CLI entrypoint loads `SpaceAGORAVisualizer` for these
commands only. When calling `run_cli` from Julia to view saved results, load both
`SpaceAGORA` and `SpaceAGORAVisualizer` first.

A visualized child uses the active project unless `SPACEAGORA_VISUALIZER_ENV`
explicitly selects another prepared project. The selected project must declare
both packages. The child loads the renderer before the example starts, so an
example can still activate the core project without losing viewer support.
The script's `PROGRAM_FILE`, arguments, and main guard retain their usual
meaning. A missing project or package declaration is rejected before launching
the example, including in `--print-only` mode. If declarations exist but the
renderer has not been installed, its preload fails before the example can
create outputs. The CLI never installs packages.

For the convenience wrappers, which start under the core project, set
`SPACEAGORA_VISUALIZER_ENV` to the prepared project directory. See
[Visualization](user/visualization.md) for the saved-payload workflow and
renderer options.

### Run telemetry verification

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

Profiles:

| Profile | Effect |
|---|---|
| `quick` | Development-scale telemetry verification |
| `full` | Full telemetry verification profile |
| `smoke` | CLI smoke alias for the quick profile |

The profile can also be set explicitly:

```text
julia --project=. src/interfaces/cli/main.jl telemetry --profile=full --output-dir=output/telemetry_full
```

Supported options:

| Option | Effect |
|---|---|
| `--output-dir=<dir>` | Writes telemetry CSV outputs under this directory |
| `--enforce=0|1` | Fails the command when verification thresholds fail |
| `--plots=0|1` | Enables or disables telemetry plot generation |
| `--scenarios=a,b` | Runs only the named manifest scenarios (default: all); unknown names are an error |
| `--print-only` | Prints the resolved launcher without running it |

Plot generation is off by default in the CLI path so the telemetry command
remains usable even when plotting assets/scripts are unavailable. Enable it
explicitly with `--plots=1`.

The telemetry CLI writes:

```text
telemetry_orbit_accuracy_summary.csv
telemetry_orbit_accuracy_errors.csv
```

### Run benchmark launchers

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_cli
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
```

Supported benchmark modes:

| Mode | Launcher |
|---|---|
| `runtime-analysis` | `benchmarks/studies/performance_runtime_analysis.jl` |
| `smart-parallel-ladder` | `benchmarks/studies/performance_smart_parallel_ladder.jl` |

Supported options:

| Option | Effect |
|---|---|
| `--output-dir=<dir>` | Sets the launcher-specific output directory environment variable |
| `--print-only` | Prints the resolved launcher without running it |

Any other arguments after the benchmark mode are passed through to the study
launcher. For example:

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis quick --output-dir=output/perf_quick
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder full --passes=3 --output-dir=output/smart_ladder_full
```

### Inspect assets

```text
julia --project=. src/interfaces/cli/main.jl assets check
julia --project=. src/interfaces/cli/main.jl assets manifest
julia --project=. src/interfaces/cli/main.jl assets setup-open
```

## Print-only mode

Every execution command supports `--print-only` so you can inspect the resolved
launcher, environment, and project without running the workload.

Examples:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --smoke --print-only
julia --project=. src/interfaces/cli/main.jl telemetry smoke --output-dir=output/telemetry_smoke --print-only
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke --print-only
```

### Named surrogate data

After [setting up the public Odyssey environment](tutorials/odyssey_surrogate.md):

```sh
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets list
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets fetch --preset=odyssey_p20_frozen_v1 --version=1.0.0
julia --project=examples/odyssey_surrogate_env src/interfaces/cli/main.jl assets check --preset=odyssey_p20_frozen_v1 --version=1.0.0
```

`list` reads the catalog without downloading. `fetch` downloads only the selected
version and verifies its bytes. The same commands take
`--preset=mars_global_upper_p20_frozen_v1 --version=1.0.0` for the global upper
atmosphere preset and `--preset=mars_global_near_surface_p20_frozen_v1 --version=1.2.0`
(or `--version=1.1.0` or `--version=1.0.0`) for the near-surface preset.
`--preset=mars_global_descent_p20_frozen_v1 --version=1.0.0` selects the composite descent preset: `fetch` installs
and verifies both of its components, `check` builds the composite from them, and it takes no `--file`, since it has
no payload of its own. `check` uses installed data offline and validates
the payload schema, coordinates or domain limits, and generation settings. `--file=<path>` selects
an explicit exact-byte copy; an invalid path never falls back to retrieval.
`fetch --offline` verifies an installed artifact without network access.
Both `--preset=<id>` and `--preset <id>` forms are accepted, likewise `version`
and `file`. An exact version is required.

The native-free wrapper loads automatically for the preset `check` command.
These commands do not build or initialize native GRAM. Normal Odyssey example
startup performs preset retrieval itself, so the commands above are optional
inspection/prefetch tools. The unqualified `assets check` reports the readiness
of repository-local asset roots and catalogs (usable, missing, empty, Git LFS
pointers, or empty payloads) and exits with status 1 when a required
entry is not usable; it confirms that payload bytes exist, not that a grid is
valid.
The printed `runtime_wrapper_revision` is the GRAMSuite revision used when the
preset was generated; the environment's own pin is in its `Project.toml`.
