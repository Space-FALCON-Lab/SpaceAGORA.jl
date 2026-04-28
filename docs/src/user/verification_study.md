# Verification Study

Use this page when you want the packaged telemetry verification workflow rather
than a one-off example run.

This page is for users validating a local setup, comparing outputs, or running
the repository-owned verification study with explicit enforcement.

Shortest successful command:

```bash
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

What to read next:

- [Assets & Modes](../assets.md)
- [CLI](../cli.md)
- [Recipes](recipes.md)

## Script entrypoint

The direct study launcher is:

```bash
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

Use `quick` when you want the shortest validation path. Keep `--enforce=true`
when threshold failures should stop the run.

## CLI entrypoint

The packaged CLI route is:

```bash
./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

Plot generation is off by default in the CLI path so the command remains usable
even when plotting dependencies or scripts are unavailable. Enable it explicitly
with `--plots=1` when needed.

## What this workflow is good for

Use the verification study when you need:

- a known repository-owned validation path
- explicit enforcement behavior
- a repeatable output directory for local inspection
