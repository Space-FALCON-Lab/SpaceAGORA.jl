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
multi-process benchmark and study workflows, plus deterministic parallel-policy
hint persistence.

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

The following are not yet treated as a full cluster product surface:

- arbitrary scenario execution on schedulers
- distributed checkpoint/restart orchestration
- automatic artifact staging across nodes

Those remain future work and should not be assumed stable unless they are
promoted into the root package API and documented here.
