# CLI

Use this page when you want the stable command-line surface instead of calling
Julia scripts directly.

This page is for users and operators who prefer packaged commands for examples,
verification studies, benchmarks, and asset inspection.

Shortest successful command:

```text
julia --project=. src/cli/main.jl assets check
```

What to read next:

- [Quickstart](user/quickstart.md)
- [Recipes](user/recipes.md)
- [Assets and Modes](assets.md)

SpaceAGORA exposes a package-owned CLI surface through:

```text
julia --project=. src/cli/main.jl
```

This Julia entrypoint is the cross-platform command path for Windows, Linux,
and macOS.

Convenience wrappers:

- Linux/macOS: `./bin/spaceagora`
- Windows: `bin\spaceagora.bat`

## Commands

### Run an example

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

Supported options:

| Option | Effect |
|---|---|
| `--example=<file>` | Example filename under `examples/`, or a direct path |
| `--output-dir=<dir>` | Sets `SPACEAGORA_CLI_OUTPUT_DIR` for the child process |
| `--smoke` | Sets the example smoke-mode environment variables |
| `--print-only` | Prints the resolved launcher without running it |

Use `--smoke` for long examples that support
`SPACEAGORA_EXAMPLE_SMOKE=1`.

### Run telemetry verification

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

Profiles:

| Profile | Effect |
|---|---|
| `quick` | Development-scale telemetry verification |
| `full` | Full telemetry verification profile |
| `smoke` | CLI smoke alias for the quick profile |

The profile can also be set explicitly:

```text
julia --project=. src/cli/main.jl telemetry --profile=full --output-dir=output/telemetry_full
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
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_cli
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
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
julia --project=. src/cli/main.jl benchmark runtime-analysis quick --output-dir=output/perf_quick
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder full --passes=3 --output-dir=output/smart_ladder_full
```

### Inspect assets

```text
julia --project=. src/cli/main.jl assets check
julia --project=. src/cli/main.jl assets manifest
julia --project=. src/cli/main.jl assets setup-open
```

## Print-only mode

Every execution command supports `--print-only` so you can inspect the resolved
launcher, environment, and project without running the workload.

Examples:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --smoke --print-only
julia --project=. src/cli/main.jl telemetry smoke --output-dir=output/telemetry_smoke --print-only
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke --print-only
```
