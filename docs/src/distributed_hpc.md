# Distributed and HPC

Use this page when you are running SpaceAGORA benchmark or study workflows with
process workers or scheduler launchers.

This page is for operators and advanced users who need reproducible multi-process
execution, worker-project setup, or scheduler examples.

Shortest successful example:

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/hpc_local_runtime_analysis
```

What to read next:

- [CLI](cli.md)
- [Parallel Execution](user/parallel_execution.md)
- [Recipes](user/recipes.md)
- [Maintainer Overview](maintainer/index.md)

SpaceAGORA's supported distributed execution surface today is centered on
multi-process benchmark and study workflows, self-bootstrapping process-backend
outer parallelism for library-level campaigns, and deterministic
parallel-policy hint persistence.

This is not a promise that every runtime path is cluster-hardened. It is the
documented path that currently has real ownership in the repository.

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

For library-level campaigns (`run_monte_carlo`, `run_constellation_ensemble`),
process-worker parallelism does not require any of the manual `addprocs` or
scheduler setup below. Passing `threads=:auto` lets the adaptive outer-route
bandit choose the `:process` route itself for workload shapes where it wins,
and auto-bootstraps a `Distributed` worker pool through
`SpaceAGORA.ParallelProcess`:

```julia
result = run_monte_carlo(1:100; threads=:auto) do seed
    run_simulation(make_config_for_seed(seed); return_solution=true)
end
```

Under the hood this uses `campaign_process_pool()` (the process-global
`ProcessPool` shared across campaigns in the session) and
`ensure_process_workers!(pool, n; warmup_fn=...)`, which spawns only the
worker shortfall, bootstraps each new worker with `SpaceAGORA`, `GRAMSuite`
(best-effort), and the default SPICE kernel set, and optionally pays a
representative call's JIT/specialization cost once per new worker (measured
at roughly 70 s cold vs. a fraction of a second warm) via `warmup_fn` instead
of inside the first real, timed dispatch. Each process worker runs with
`--threads=1`, so it does not contend with the coordinator's thread pool.
`shutdown_process_pool!(pool)` tears the pool down; campaign code otherwise
leaves it warm across calls by design.

See [Parallel Execution: Process-backend outer parallelism](user/parallel_execution.md#process-backend-outer-parallelism)
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

The Slurm launcher shows:

- project selection for spawned workers
- explicit machine labeling
- explicit persistent-hint path selection
- a per-job writable depot path

## Current scope

The current first-class distributed support is:

- process-worker benchmark and study execution
- deterministic persistent-hint state for parallel-policy tuning
- self-bootstrapping process-backend outer parallelism for `run_monte_carlo`
  and `run_constellation_ensemble` campaigns (`SpaceAGORA.ParallelProcess`)

The following are not yet treated as a full cluster product surface:

- arbitrary scenario execution on schedulers
- distributed checkpoint/restart orchestration
- automatic artifact staging across nodes

Those remain future work and should not be assumed stable unless they are
promoted into the root package API and documented here.
