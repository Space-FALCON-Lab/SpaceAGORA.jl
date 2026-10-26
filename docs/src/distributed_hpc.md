# Distributed and HPC

SpaceAGORA supports two ways of using process workers: the benchmark and study
launchers, which take explicit worker settings, and library campaigns, whose
planner can choose a process backend. A custom campaign function needs extra
preparation before it can use process workers, described below. Arbitrary
scheduler-driven simulations, distributed checkpoint orchestration and
automatic staging of files across nodes are not supported.

From a prepared repository environment, preview a study launch from the
repository root:

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/hpc_local_runtime_analysis --print-only
```

`--print-only` shows the command without running it. The CLI starts the study
in the examples environment, so before running it for real, prepare that
environment with `julia --startup-file=no examples/setup.jl`, along with the
study's assets and worker environment, then drop `--print-only`. The
[CLI](cli.md) page lists the command options, and
[Parallel Execution](user/parallel_execution.md) covers campaign routing and
worker preparation.

## Process workers for benchmark launchers

The benchmark launchers start their process workers in
`benchmarks/studies/performance_runtime_analysis/case_catalog/env_workers_and_routes.jl`.
Three environment variables control them:

| Variable | Meaning |
|---|---|
| `SPACEAGORA_PERF_PROCS` | how many process workers a benchmark or study launcher starts |
| `SPACEAGORA_PERF_WORKER_PROJECT` | the Julia project the workers are started with |
| `SPACEAGORA_PERF_MACHINE_LABEL` | the machine label that persistent hint signatures include |

Use `.`, the root project, as the worker project unless you have a reason to
start the workers in another environment.

## Process backend for campaigns

In a library campaign, the parallel flag lets the planner choose among the
eligible routes. `run_constellation_ensemble` works out what it needs to know
about the workload from its configuration. A plain `run_monte_carlo` callback
tells the planner much less, so without explicit `route_features`,
`parallel=true` keeps the work inside the coordinator process. A process route
also needs workload metadata and workers that already have the callback's
definitions, dependencies and assets: starting the pool does not copy functions
or state from `Main` to the workers.

For example, once your study defines `make_config_for_seed`, this fragment asks
for parallel execution but does not by itself ask for process workers:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; parallel=true) do seed
    run_simulation(make_config_for_seed(seed); return_solution=true)
end
```

[Parallel Execution](user/parallel_execution.md) explains the difference
between in-process execution and the process route. The benchmark worker
settings above don't replace that preparation.

When a campaign does take the process backend, it uses
`campaign_process_pool()`, the `ProcessPool` shared by every campaign in the
session, and `ensure_process_workers!(pool, n; warmup_fn=...)`. That call starts
only the missing workers and loads each new one with `SpaceAGORA`, with
`GRAMSuite` when it can, and with the default SPICE kernels. If you pass
`warmup_fn`, it runs a representative call on each new worker before the first
dispatch, so that compilation and initialization happen up front; how much that
covers depends on the workload, the argument types, the Julia environment and
the caches, and later calls can still compile. Each worker runs with
`--threads=1`, so it does not compete with the coordinator's threads. The pool
stays up between calls until `shutdown_process_pool!(pool)` shuts it down.

[Process-backend outer parallelism](user/parallel_execution.md#Process-backend-outer-parallelism)
walks through this route. It is separate from the `addprocs`-based benchmark
route above: it needs neither `SPACEAGORA_PERF_PROCS` nor a worker project,
because the pool starts its workers in the coordinator's own active project.

Two settings decide how outer and inner parallel work share a machine's threads
once a route is chosen:

- `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER` (default `0`). The batched
  spherical-harmonics gravity route runs once per right-hand-side evaluation,
  not once per sample, so by default it runs serially whenever outer
  parallelism (`:threads` or `:process`) is active, instead of starting a nested
  batch that would compete with the outer split. Turn it on only after you have
  measured that splitting the thread budget helps your workload.
- `SPACEAGORA_DENSITY_CALLBACK_PARALLEL` (`off`, `auto` or `on`; default `auto`).
  Native GRAM is serialized behind a process-wide lock whatever the thread
  count, while the GRAM surrogate needs no lock, so under `auto` the native
  model needs a larger thread budget than the surrogate before the density
  callback uses threads.
  [Parallel Execution](user/parallel_execution.md) gives the full set of
  controls.

## Reproducible hint state

The parallel policy can keep persistent hints between runs. These variables
control it:

- `SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS`
- `SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST`
- `SPACEAGORA_PARALLEL_POLICY_STATE_RESET`
- `SPACEAGORA_PARALLEL_POLICY_STATE_PATH`
- `SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION`
- `SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES`

To keep route and hint decisions the same across repeated benchmark runs:

1. Fix the route and policy controls for the run, such as
   `SPACEAGORA_PERF_PARALLEL_BACKEND`, `SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE`
   and `SPACEAGORA_PARALLEL_POLICY_ADAPTIVE`.
2. Set a fixed `SPACEAGORA_PERF_MACHINE_LABEL`.
3. Set `SPACEAGORA_PARALLEL_POLICY_STATE_PATH`, so that repeated jobs reuse the
   same hint file on purpose rather than through the default path.
4. Use `SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` only for cold-start studies.

## Shared depot

For multi-node or multi-job runs, use a project that is already instantiated
(`--project=.`), a shared depot prepared ahead of time and mostly read, and,
where the cluster allows it, a writable depot overlay local to each job. Don't
let many concurrent jobs install packages or artifacts into the same depot for
the first time: run `Pkg.instantiate()` as a provisioning step, not inside the
benchmark job, and point `SPACEAGORA_PERF_WORKER_PROJECT` at a project that each
worker can already reach.

## Example launchers

The repository has two launcher scripts:

- `scripts/hpc/local_process_runtime_analysis.sh`, the smallest setup for one
  machine with process workers;
- `scripts/hpc/slurm_process_runtime_analysis.sh`, a starting point for a
  site-specific Slurm job script. It sets the worker project, a machine label,
  an explicit hint-file path and a writable depot path for the job.

Both are Unix shell scripts. On Windows, set the same environment variables in
PowerShell or `cmd.exe` and run the Julia CLI entry point directly.

The Slurm script's default `JULIA_DEPOT_PATH` is a new directory for the job,
so populate it before the job runs, or change the script to include your
prepared shared depot. An empty writable directory does not make the project's
dependencies available.

## Current scope

Supported distributed use is:

- process-worker benchmark and study runs;
- reproducible persistent-hint state for tuning the parallel policy;
- the self-starting process backend (`SpaceAGORA.ParallelProcess`) for
  `run_constellation_ensemble` campaigns, and for `run_monte_carlo` campaigns
  given `route_features` and prepared workers.

Running arbitrary scenarios through a scheduler, distributed checkpoint and
restart, and automatic staging of artifacts across nodes are not supported.
Don't rely on them unless they become part of the root package API and this
page documents them.
