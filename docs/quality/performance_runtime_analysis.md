# SpaceAGORA runtime analysis

The runtime-analysis study measures how long representative simulations take
to copy their configuration and to solve, under the parallel backends and
profiles below. Its launcher is `benchmarks/studies/performance_runtime_analysis.jl`,
and the CLI's `benchmark runtime-analysis` command starts the same launcher
(see [CLI](../src/cli.md)).

## Run

The study runs in the examples environment. The profile is `quick` (the
default), `full`, or `smoke`, another name for `quick`:

```bash
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl full
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl --profile=quick --outdir=output/performance_quick
```

`--outdir` overrides the output directory (default `output/performance/`); the
variables `SPACEAGORA_PERF_PROFILE` and `SPACEAGORA_PERF_OUTDIR` set the same two
things.

The study times `deepcopy(case.args_template)` on its own and then solves with
`isolate_state=false`, so the copy and setup figures it reports are the cost an
expert caller avoids by turning isolation off. Keep `isolate_state=true`, the
default, for correctness. Treat `false` as a throughput setting only when the
configuration instance is used once and is not shared with concurrent or
state-changing runs.

## Parallel backend

`SPACEAGORA_PERF_PARALLEL_BACKEND` selects `threads`, `process`, `none` or `auto`
(the default), and `SPACEAGORA_PERF_PROCS` the number of process workers:

```bash
SPACEAGORA_PERF_PARALLEL_BACKEND=threads julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
SPACEAGORA_PERF_PARALLEL_BACKEND=process SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
JULIA_NUM_THREADS=4 SPACEAGORA_PERF_PARALLEL_BACKEND=auto SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
```

With `threads` and a single Julia thread, the study warns and runs serially.
Any other value is an error.

## Layout

The launcher loads `benchmarks/studies/performance_runtime_analysis/main.jl`,
which brings in four parts, each with a folder of the same name for its pieces:

- `case_catalog.jl` (`case_catalog/`): the cases, profiles and worker and route setup;
- `measurement.jl` (`measurement/`): warm-up and batches, Monte Carlo seeds,
  entry durations, per-orbit timing, route partitioning and the rollout solve;
- `reporting.jl` (`reporting/`): summaries, report writing, rollout gates and plots;
- `cli.jl`: the command-line options.

The parallel layers the study exercises live in these files:

- the outer route: `src/orchestration/parallel/routing/parallel_profiles.jl`, and
  the launcher itself;
- the density, thermal and control callbacks:
  `src/orchestration/simulation/callbacks/density_callbacks.jl`,
  `thermal_callbacks.jl` and `control_callbacks.jl`;
- the multibody and link kernels:
  `src/models/dynamics/effectors/aerodynamic_wrench_models.jl` and
  `src/models/dynamics/effectors/perturbations.jl`;
- the solve: `src/orchestration/simulation/engine/execution.jl`.

## Outputs

The study writes, under `output/performance/` by default:

1. `runtime_raw_<profile>_<timestamp>.csv`
2. `runtime_summary_<profile>_<timestamp>.csv`
3. `runtime_per_orbit_raw_<profile>_<timestamp>.csv`
4. `runtime_per_orbit_summary_<profile>_<timestamp>.csv`
5. `runtime_report_<profile>_<timestamp>.md`
6. `runtime_plot_*_<profile>_<timestamp>.png`

These files belong in ignored `output/` paths. If CI needs them, publish them
as workflow artifacts; don't commit machine-specific reports to the repository.

## Retry records and accounting

The runtime, per-orbit and entry-duration raw tables keep one row for every
measured attempt that returned, failures before a successful retry included.
The scenario, seed, repeat, sweep and reference-or-measured fields identify the
request, and `attempt` numbers its attempts in order. Collection stops at the
first success or at the configured attempt limit, and `is_terminal_attempt`
marks that final outcome. Adaptive routing feedback uses only terminal
outcomes; rows from older files without the marker keep their earlier meaning.
`measure_montecarlo_seed`, for a single seed, returns the terminal row and the
error message.

In the summaries, `samples_total`, `samples_failed`, the success and failure
rates and the Success/Total tables count attempts, while the timing statistics
use successful rows only. In a complete, newly collected sequence, rows with
`attempt == 1` count the requested runs and the other rows count retries, so
the attempt-level rates are not completion rates of requested runs. Older files
that kept only terminal rows cannot recover the discarded attempts or their
costs.

`total_time_all_attempts_s` sums the recorded copy-plus-solve costs, and
`penalized_expected_wall_time_s` divides that sum by the number of requested
runs. Neither is the wall time of a whole campaign: explicit garbage
collection, warm-ups, worker startup and later processing are not measured. An
interrupt or a strict error can stop the study before a row is returned or the
outputs are written, and nothing recovers measurements that were never
returned or written.

## Reporting ownership and isolated tests

`benchmarks/studies/performance_runtime_analysis/reporting/runtime_reports.jl`
owns the profile record, the result summaries, the text formatting, the
trajectory comparison and the rollout decisions, and can be imported on its
own. It needs CSV, DataFrames and standard libraries, but loading it does not
load SpaceAGORA, GRAMSuite, SPICE or Plots, start workers, activate a project
or change environment settings; settings are read when its functions are
called, not at import.

The simulation entry point imports this module before it builds the case
catalog. `reporting.jl` gives the report writer the current machine, solver and
guardrail context, and gives the rollout evaluators their warm-up and
simulation callbacks. The native simulation runner is in
`measurement/rollout_solution.jl`, and plotting is a separate component. A
writer used on its own must supply the report context explicitly; provenance
recorded in a row takes precedence over the fallback machine fields.

The reporting tests import the real module. They share a synthetic row builder,
`test/support/runtime_report_fixtures.jl`, and keep their expected values,
failure controls and CSV and Markdown assertions in the individual tests. A
fresh-process test checks that importing the module leaves the environment,
project, load path, working directory and worker list unchanged. These tests
check reporting contracts; they say nothing about simulation accuracy or speed.

### Shared benchmark support

`benchmarks/support/artifacts.jl` finds benchmark artifacts. The static
comparison wrapper and the smart ladder both use its newest-modification-time
lookup: a missing required artifact is an error, and an optional lookup returns
`nothing` only when the output directory exists but holds no matching artifact.

The entry-duration and per-orbit collectors share `measurement/routing.jl`,
which partitions the work into ordered routes and falls back to serial
execution when a thread request is unsafe. Each collector keeps its own
execution, reference runs, row and log assembly and environment restoration.
Importing these support modules starts no workers and no simulation. Tests
cover the routing decisions, artifact failures and the retry and output-order
contracts. The runtime-Any gate tests build their fixtures with
`test/support/runtime_any_gate_fixtures.jl` and keep their own negative
controls in the individual test files.

## Standalone parallelization performance study

`benchmarks/studies/parallelization_performance.jl` is a separate study from
the runtime-analysis and smart-ladder harnesses, with its own case catalog,
mode matrix, subprocess runner, trajectory parity checks, CSV summaries, report
generation and paper protocol script. Its profiles are `smoke` (the default),
`full` and `test`:

```bash
julia --project=examples benchmarks/studies/parallelization_performance.jl smoke
julia --project=examples benchmarks/studies/parallelization_performance.jl full --outdir=output/performance/parallelization_performance_full
bash benchmarks/studies/parallelization_performance/protocol.sh
```

The last command is the paper protocol. The study writes:

1. `parallelization_performance_raw_<profile>_<timestamp>.csv`
2. `parallelization_performance_summary_<profile>_<timestamp>.csv`
3. `parallelization_trajectory_parity_<profile>_<timestamp>.csv`
4. `parallelization_hardware_<profile>_<timestamp>.csv`
5. `parallelization_performance_report_<profile>_<timestamp>.md`

A speedup claim should use only rows whose matching serial baseline completed
and whose trajectory parity probe passed, where the scenario configures one.
