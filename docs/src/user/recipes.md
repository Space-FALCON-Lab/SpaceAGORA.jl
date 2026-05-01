# Recipes

Use this page when you want copy-paste commands for common tasks.

This page is for users who already know the task they want to perform and do
not need a longer explanation first.

Shortest successful command:

```text
julia --project=. src/cli/main.jl assets check
```

What to read next:

- [Quickstart](quickstart.md)
- [CLI](../cli.md)
- [Verification Study](verification_study.md)

## Run the first no-GRAM example

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

## Run an example through the CLI

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

## Inspect local assets

```text
julia --project=. src/cli/main.jl assets check
julia --project=. src/cli/main.jl assets manifest
```

## Run the telemetry verification study

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

## Run the telemetry verification CLI path

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

## Build the docs locally

```text
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
julia --project=docs docs/make.jl
```

## Regenerate common local artifacts

```text
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```
