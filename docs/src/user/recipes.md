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
- [Examples Catalog](examples_catalog.md)
- [Studies and Benchmarks](studies_benchmarks.md)
- [Verification Study](verification_study.md)
- [Simulation Outputs](outputs.md)

## Run the first no-GRAM example

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

## Run an example through the CLI

```text
julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

## Run an aerobraking example in smoke mode

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

## Run a Monte Carlo example script with threads

```text
julia --project=. --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

Inside the script:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args)
end
```

## Run an RPO planner-comparison smoke case

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=. examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

## Run the robot-arm and Cloth dynamics smoke batch

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=. examples/Robot_Arm_Planner_Cloth_Demo.jl
```

## Inspect local assets

```text
julia --project=. src/cli/main.jl assets check
julia --project=. src/cli/main.jl assets manifest
```

## Run the telemetry verification study

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
# one scenario only (also honoured as SPACEAGORA_TELEMETRY_SCENARIOS=odyssey)
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true --scenarios=odyssey
```

## Run the telemetry verification CLI path

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

## Run benchmark smoke launchers

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_smoke
```

## Inspect a benchmark command before running it

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis quick --output-dir=output/perf_quick --print-only
```

## Build the docs locally

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. docs/make.jl
```

## Regenerate common local artifacts

```text
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```
