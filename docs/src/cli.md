# CLI

Use this page when you want the stable command-line surface instead of calling
Julia scripts directly.

This page is for users and operators who prefer packaged commands for examples,
verification studies, benchmarks, and asset inspection.

Shortest successful command:

```bash
./bin/spaceagora assets check
```

What to read next:

- [Quickstart](user/quickstart.md)
- [Recipes](user/recipes.md)
- [Assets and Modes](assets.md)

SpaceAGORA exposes a package-owned CLI surface through:

```bash
./bin/spaceagora
```

The wrapper keeps execution on the committed `/.AGORA` environment by default
and forwards into `SpaceAGORA.run_cli(...)`.

## Commands

### Run an example

```bash
./bin/spaceagora run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
```

Use `--smoke` to force the example smoke configuration.

### Run telemetry verification

```bash
./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

Plot generation is off by default in the CLI path so the telemetry command
remains usable even when plotting assets/scripts are unavailable. Enable it
explicitly with `--plots=1`.

### Run benchmark launchers

```bash
./bin/spaceagora benchmark runtime-analysis smoke --output-dir=output/perf_cli
./bin/spaceagora benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
```

### Inspect assets

```bash
./bin/spaceagora assets check
./bin/spaceagora assets manifest
./bin/spaceagora assets setup-open
```

## Print-only mode

Every execution command supports `--print-only` so you can inspect the resolved
launcher, environment, and project without running the workload.
