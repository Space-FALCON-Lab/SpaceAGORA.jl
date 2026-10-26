# Distributed and HPC

This guide is for users preparing process workers or scheduler jobs for
SpaceAGORA studies. It covers two routes: benchmark launchers with explicit
worker settings, and library campaigns whose planner can select a process
backend. Custom campaign functions need additional preparation, described below.

Start with a prepared repository environment and preview the study command from
the repository root:

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/hpc_local_runtime_analysis --print-only
```

The preview shows the launch command without running the study. The CLI starts
the study in the examples environment, so before executing it, prepare that
environment with `julia --startup-file=no examples/setup.jl`, and the selected
study's assets and worker environment; remove `--print-only` when those
prerequisites are ready. See the [CLI](cli.md) for
command options and [Parallel Execution](user/parallel_execution.md) for campaign
routing and worker preparation.

The maintained distributed workflows are benchmark and study launches and
process-backed campaigns. Arbitrary scheduler-driven simulations, distributed
checkpoint orchestration and automatic staging across nodes are outside the
scope described here.

## Supported process-worker path

The canonical process-worker bootstrap logic lives in:

- `benchmarks/studies/performance_runtime_analysis/case_catalog/env_workers_and_routes.jl`

The key environment knobs are:

- `SPACEAGORA_PERF_PROCS`
  - target process-worker count for benchmark/study launchers
- `SPACEAGORA_PERF_WORKER_PROJECT`
  - Julia project path used when spawning process workers
- `SPACEAGORA_PERF_MACHINE_LABEL`
  - machine label folded into persistent hint signatures

Use `.` (the root project) as the worker project unless you have a deliberate
reason to point workers at a different environment.

## Campaign-level process backend

For library-level campaigns, the parallel flag lets the planner select among
eligible routes. `run_constellation_ensemble` derives workload information from
its configuration. A generic `run_monte_carlo` callback supplies less information:
without explicit `route_features`, `parallel=true` keeps its candidates within
the coordinator process. A process route additionally needs workload metadata
and workers prepared with the callback's definitions, dependencies and assets.
Starting the pool does not make arbitrary functions or state in `Main` available
on workers.

For example, after your study defines `make_config_for_seed`, this fragment
requests parallel execution but does not by itself request process workers:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; parallel=true) do seed
    run_simulation(make_config_for_seed(seed); return_solution=true)
end
```

See [Parallel Execution](user/parallel_execution.md) for the distinction between
in-process execution and the advanced process route. The benchmark worker
settings described above do not replace that campaign preparation.

When a campaign selects the process backend, it uses `campaign_process_pool()` (the process-global
`ProcessPool` shared across campaigns in the session) and
`ensure_process_workers!(pool, n; warmup_fn=...)`, which spawns only the
worker shortfall, bootstraps each new worker with `SpaceAGORA`, `GRAMSuite`
(best-effort), and the default SPICE kernel set, and optionally pays a
representative call's compilation and initialization cost via `warmup_fn`
before the first timed dispatch. Warmup is attempted on each new worker; its
cost and coverage depend on the workload, argument types, Julia environment
and cache state. Later calls can still require additional compilation. Each
process worker runs with
`--threads=1`, so it does not contend with the coordinator's thread pool.
`shutdown_process_pool!(pool)` tears the pool down; campaign code otherwise
leaves it warm across calls by design.

See [Parallel Execution: Process-backend outer parallelism](user/parallel_execution.md#Process-backend-outer-parallelism)
for the full walkthrough. This path is entirely separate from the
scheduler/`addprocs`-based benchmark path documented above — it needs no
`SPACEAGORA_PERF_PROCS` or worker-project configuration, since the pool
bootstraps itself against the coordinator's own active project.

Two related tuning knobs affect how outer- and inner-parallel work share a
machine's threads once a route is chosen:

- `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER` (default `0`) — the
  spherical-harmonics gravity SIMD batch route fires once per RHS/ODE step
  rather than once per sample, so by default it runs serially whenever outer
  parallelism (`:threads` or `:process`) is active, rather than spawning a
  nested worker batch that would contend with the outer split. Only enable
  this after measuring that splitting the thread budget helps your workload.
- `SPACEAGORA_DENSITY_CALLBACK_PARALLEL` (`off`/`auto`/`on`) — the density
  callback's inner-parallelism thread floor is lower for the lock-free GRAM
  surrogate than for native GRAM (which serializes through a process-wide
  lock regardless of thread count), so `auto` picks a different effective
  floor per density-model family. See
  [Parallel Execution](user/parallel_execution.md) for the full control
  matrix.

## Deterministic hint state

The parallel-policy persistent hint layer is controlled by:

- `SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS`
- `SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST`
- `SPACEAGORA_PARALLEL_POLICY_STATE_RESET`
- `SPACEAGORA_PARALLEL_POLICY_STATE_PATH`
- `SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION`
- `SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES`

Recommended practice:

1. Set stable route and policy controls for the run, such as
   `SPACEAGORA_PERF_PARALLEL_BACKEND`,
   `SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE`, and
   `SPACEAGORA_PARALLEL_POLICY_ADAPTIVE`.
2. Set a stable `SPACEAGORA_PERF_MACHINE_LABEL`.
3. Set an explicit `SPACEAGORA_PARALLEL_POLICY_STATE_PATH` so repeated jobs
   reuse the same hint file intentionally instead of relying on the default
   path.
4. Use `SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` only for cold-start studies.

That keeps route and hint decisions reproducible across repeated benchmark runs.

## Shared depot guidance

For multi-node or multi-job execution, prefer:

1. a pre-instantiated project (`--project=.`)
2. a read-mostly shared depot that is prepared ahead of time
3. a job-local writable depot overlay when cluster policy allows it

Practical guidance:

- do not rely on many concurrent jobs mutating the same depot for first-time
  package/artifact installation
- run `Pkg.instantiate()` as a provisioning step, not inside the benchmark job
- keep `SPACEAGORA_PERF_WORKER_PROJECT` pointed at a project already available
  on each worker

## Example launchers

Repository examples:

- `scripts/hpc/local_process_runtime_analysis.sh`
- `scripts/hpc/slurm_process_runtime_analysis.sh`

The local launcher is the minimal path for one machine with process workers.
Those helper scripts are Unix-oriented examples; on Windows, set the same
environment variables in PowerShell or `cmd.exe` and run the Julia CLI
entrypoint directly.

The Slurm launcher is a starting point for a site-specific job script. Its
default `JULIA_DEPOT_PATH` is a job-local path, so it must be populated before
execution or adapted to include the prepared shared depot. Merely creating an
empty writable directory does not make the project dependencies available.

The Slurm launcher shows:

- project selection for spawned workers
- explicit machine labeling
- explicit persistent-hint path selection
- a per-job writable depot path

## Current scope

The current first-class distributed support is:

- process-worker benchmark and study execution
- deterministic persistent-hint state for parallel-policy tuning
- self-bootstrapping process-backend outer parallelism
  (`SpaceAGORA.ParallelProcess`) for `run_constellation_ensemble` campaigns, and
  for `run_monte_carlo` campaigns given `route_features` and prepared workers

The following are not yet treated as a full cluster product surface:

- arbitrary scenario execution on schedulers
- distributed checkpoint/restart orchestration
- automatic artifact staging across nodes

Those remain future work and should not be assumed stable unless they are
promoted into the root package API and documented here.
