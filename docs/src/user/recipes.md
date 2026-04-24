# Recipes

Use this page when you want copy-paste commands for common tasks.

This page is for users who already know the task they want to perform and do
not need a longer explanation first.

Shortest successful command:

```bash
./bin/spaceagora assets check
```

What to read next:

- [Quickstart](quickstart.md)
- [CLI](../cli.md)
- [Verification Study](verification_study.md)

## Run the first no-GRAM example

```bash
julia --project=.AGORA examples/AGORA_Earth_NoGRAM.jl
```

## Run an example through the CLI

```bash
./bin/spaceagora run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
```

## Inspect local assets

```bash
./bin/spaceagora assets check
./bin/spaceagora assets manifest
```

## Run the telemetry verification study

```bash
julia --project=.AGORA benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

## Run the telemetry verification CLI path

```bash
./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

## Build the docs locally

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Regenerate common local artifacts

```bash
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl docs
```
