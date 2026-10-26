# Studies and Benchmarks

Use this page when you need a repeatable study or benchmark, not just one
example run.

This page is for users generating verification reports, runtime measurements,
GRAM comparison outputs, or paper-study artifacts.

Inspect the benchmark command without starting a workload:

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke --print-only
```

What to read next:

- [CLI](../cli.md)
- [Parallel Execution](parallel_execution.md)
- [Distributed and HPC](../distributed_hpc.md)
- [Simulation Outputs](outputs.md)

## Start here

After checking the prerequisites below, run a reduced benchmark workload:

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke
```

Benchmark `smoke` profiles execute workloads, including simulations. Executing
the `runtime-analysis` or `smart-parallel-ladder` launcher needs the
`data/GRAMSuite.jl` submodule ([GRAMSuite Setup](gramsuite_setup.md)) because
both import GRAMSuite before parsing study options. A prepared caller package
takes precedence; the vendored project is an append-only discovery fallback.
The maintained runtime analysis uses package-loaded SpaceAGORA models and the
canonical GRAM extension. Its native-backed results use a different numerical
route from historical raw-loaded measurements; see the
[migration evidence](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/137b55d43066054ab9534c6d9123ce9b5dffd20d/docs/quality/runtime_benchmark_loading_migration.md). GRAM-backed
cases also need the native GRAM library built. The CLI runs its child in the
examples environment; prepare it with `julia --startup-file=no examples/setup.jl`
before execution.

Paper and mission studies use the separate repositories listed below. Their
source pins, environments and data prerequisites belong to those projects.

If you need telemetry verification outputs, run:

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

If you need parallel-policy comparison outputs, run:

```text
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_smoke
```

Use `--print-only` on any CLI-backed study to inspect the resolved Julia command
before running it.

## Choose by output

### Telemetry verification report

Use this when the goal is orbit-event or time-aligned telemetry validation.

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

The command writes:

```text
telemetry_orbit_accuracy_summary.csv
telemetry_orbit_accuracy_errors.csv
```

Add `--plots=1` when plot artifacts are required.

### Runtime analysis report

Use this when the goal is timing, scaling, hardware metadata, or performance
regression artifacts.

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis quick --output-dir=output/performance
```

Use `smoke` first on a new machine. Use `full` only when you need full study
artifacts.

### Smart parallel ladder

Use this when the goal is comparing serial, threaded, process, and adaptive
parallel-policy rungs.

```text
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
```

## Profiles

| Profile | Use for |
|---|---|
| `smoke` | Reduced workload after its prerequisites are ready |
| `quick` | Development and regression artifacts |
| `full` | Full study or paper-grade artifacts |

## Research repositories

Use these repositories for the research campaigns moved out of core. Each records
its original source, pins its SpaceAGORA dependency and keeps its own environment,
outputs and study tests. Follow that repository's setup instructions before
running a script; core's examples environment does not prepare a research project.

| Research project | What it contains |
|---|---|
| Atmosphere studies (SpaceAGORAAtmosphereStudies) | GRAM fidelity comparisons, aerobraking perturbation Monte Carlo and historical entry-density investigations |
| Parallelization studies (SpaceAGORAParallelizationStudies) | Paper protocols, constellation and Monte Carlo scaling investigations, cross-machine comparisons and paper-result tooling |
| Mission reconstructions (SpaceAGORAMissionReconstructions) | CYGNSS reconstruction, telemetry fitting, mission-specific scenario sweeps and dispersed propagation investigations |
| RPO studies (SpaceAGORARPOStudies) | The final iSpaRo HyPR reconstruction, its saved inputs and study-specific checks |

A completed run is not by itself evidence of physical agreement or reproduction
of a publication. Read each project's interpretation limits. Some studies need
licensed native libraries, existing spacecraft assets or separately obtained
telemetry. The study repositories are private. Lab members retain access; outside
collaborators and students need an invitation. Dataset access and licensing
requirements still apply independently.

The former `SpaceAGORACalibration.jl/` campaign package and its `bin/` drivers
remain retired. Use the mission repository for historical telemetry fitting, or
[Verification Study](verification_study.md) for the maintained core workflow that
compares a simulation with observations. Automatic performance tuning for
parallel execution remains available in core.

## What stays in core

The CLI entry points above remain here, together with small teaching examples,
numerical regression fixtures and reusable verification and benchmark harnesses.
Their implementations include `performance_runtime_analysis/`,
`parallelization_performance/` and the smart parallel ladder under
`benchmarks/studies/`. A folder named `studies` does not make every file a paper
campaign: these retained tools support core development and its CLI.

Core tests retain the generic planner, actuator, telemetry-parser and parallel
execution assertions. Tests tied to a particular publication, saved paper input
or mission record move with that research project. The replacement projects run
their own checks; core does not silently fetch their studies to run its tests.

Shared installation, calibration and verification commands remain in `scripts/`.
Development probes live in `scripts/dev/`. Core's remote launcher still prepares
and runs core jobs; the paper workload and its remote recipe belong to the
parallelization repository.

## Record a result

Set an explicit output directory and retain the source revision, environment,
configuration, host, logs and raw outputs together. Keep generated plots, profiler
captures and result tables outside the source checkout. A historical result does
not validate current code.

The paper-routing analysis and archive helper now live beside their runner in the
parallelization repository. Its `CASES.md` explains which tools consume raw
per-repeat records and which consume aggregated successful runs. Use that
project's archive destination rather than writing results back into core.
