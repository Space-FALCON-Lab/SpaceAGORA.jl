# Native-GRAM worker memory model

Before the outer router offers the process route,
`src/orchestration/parallel/routing/machine_topology.jl` estimates how much
memory each process worker will hold, and `effective_process_workers` and
`mixed_local_slots` (in `outer_route_selection.jl`) cap the route by what fits.
The estimate has a worker floor and, for native-GRAM point-density workloads, a
per-spacecraft term. It is a routing estimate, not a guarantee that a given
mission fits in memory. The measurements behind the current constants, and
earlier fits, are kept in the lab's research records
(see [Record a result](../src/user/studies_benchmarks.md#record-a-result)).

## Current model

| Term | Value | Where |
|---|---|---|
| Worker floor | 2 GiB, or this process's resident set if larger | `_WORKER_MEMORY_FLOOR_BYTES` |
| Per spacecraft, direct native point density | 2 MiB | `_GRAM_SAT_MEMORY_BYTES_BY_PATH.point` |
| Per spacecraft, density frozen per step | 2 MiB | `_GRAM_SAT_MEMORY_BYTES_BY_PATH.freeze_per_step` |
| Per spacecraft, look-ahead cache | 2 MiB | `_GRAM_SAT_MEMORY_BYTES_BY_PATH.lookahead` |
| Per spacecraft, offline surrogate | 0 | `_GRAM_SAT_MEMORY_BYTES_BY_PATH.surrogate` |

The measurements did not separate the three native paths, so each is charged
the largest measured slope with a safety margin. They covered short arcs with
the study's save settings and solutions not returned; longer missions, denser
saving or samples that keep their solutions may cost more. On a tight memory
budget, set `SPACEAGORA_GRAM_SAT_MEMORY_MB` (the per-spacecraft term for every
path at once) or `SPACEAGORA_PERF_WORKER_MEMORY_GB` (the worker's base
estimate) from a measured worker footprint.

`test/unit/parallel/gram_memory_model_tests.jl` checks the relations between
these terms and the routing arithmetic they feed, on fabricated machine facts.

The per-path table is kept even though its entries are equal, so that a future
per-path measurement changes one entry. Today it cannot change routing on its
own: `OuterRouteFeatures.density_family` does not distinguish the point, frozen
and look-ahead modes (they are selected by environment variables), and both call
sites of `native_gram_worker_extra_bytes` use the default `path=:point`. A fit
that separated the paths would also need a way to carry the chosen path to those
call sites.

## Measurement method

`benchmarks/studies/gram_memory_footprint/run_footprint.jl` launches a separate
one-thread Julia process (`footprint_worker.jl`) for each density path and
constellation size. That matches the shape of a process-pool worker and keeps
one point's heap out of the next point's baseline. A one-spacecraft warm-up pays
the first solver and native allocations before the N-spacecraft configuration
is measured.

The worker records peak RSS, live RSS before and after a full collection, and
Linux `VmRSS` when available. `fit_footprint.jl` fits an intercept and a
per-spacecraft slope to the peak readings and does not subtract the baseline
RSS, because the kernel may not have released pages after the collection when
the baseline is read. Include repeated solves (`GMF_SOLVE_REPEATS`), because
pool workers serve more than one sample.

| Path | Configuration |
|---|---|
| `none` | no atmosphere |
| `gram_point` | direct native point density, no freeze or vacuum cache |
| `gram_freeze_step` | density frozen per accepted step |
| `gram_lookahead` | vacuum-predicted look-ahead cache, with the horizon and deviation limits extended past the mission |
| `gram_surrogate` | offline surrogate table |

The study flies an Earth LEO constellation with the same geometry, force stack,
tolerances and save settings as the parallelization paper scenarios (kept in the
SpaceAGORAParallelizationStudies repository), and leaves state isolation at the
production default. Record any setting you change with the new fit.

```bash
julia --project=. benchmarks/studies/gram_memory_footprint/run_footprint.jl
julia --project=. benchmarks/studies/gram_memory_footprint/fit_footprint.jl benchmarks/studies/gram_memory_footprint/results/gram_memory_footprint_<host>.csv
julia --project=. --threads=8 benchmarks/studies/gram_memory_footprint/route_fit_table.jl
```

`route_fit_table.jl` prints the per-worker charge, the affordable process
workers and the mixed-dispatch local slots the model produces for machine
profiles supplied through `SPACEAGORA_CORE_BUDGET`,
`SPACEAGORA_MEMORY_BUDGET_GB`, `SPACEAGORA_MEMORY_AVAILABLE_GB` and
`SPACEAGORA_PERF_WORKER_MEMORY_GB`, so it answers for the profile rather than
the host it runs on.

Study controls are `GMF_PATHS`, `GMF_SIZES`, `GMF_MISSION_S`, `GMF_GRAVITY`,
`GMF_ISOLATE`, `GMF_REPEATS`, `GMF_SOLVE_REPEATS`, `GMF_TIMEOUT_S`,
`GMF_MIN_FREE_GB` and `GMF_SUFFIX`; `run_footprint.jl` documents each one and
its default. Results go to the study's ignored `results/` directory. Keep the
raw rows and the host and revision metadata with the fit in the private research
record, not in this repository.

## Limits

The fit does not price a density-service pool or isolated native pools inside
one process, and it does not model inner-thread width. Long missions, denser
output and other force stacks need their own measurement. Surrogate queries must
use a table built for the scenario epoch, and a surrogate run says nothing about
native worker memory. Timings taken on a shared host are diagnostic, not
portable performance figures.
