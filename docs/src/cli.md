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

Use `--smoke` to force the example smoke configuration.

### Run telemetry verification

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

Plot generation is off by default in the CLI path so the telemetry command
remains usable even when plotting assets/scripts are unavailable. Enable it
explicitly with `--plots=1`.

### Run benchmark launchers

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_cli
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
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
