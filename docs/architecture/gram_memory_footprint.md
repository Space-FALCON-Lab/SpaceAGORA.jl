# Native-GRAM worker memory model

`src/orchestration/parallel/routing/machine_topology.jl` estimates each process
worker's memory footprint before `effective_process_workers` and
`mixed_local_slots` admit a process route. The model includes a worker floor
and a per-spacecraft term for native-GRAM point-density workloads. These are
routing estimates, not a guarantee that an arbitrary mission fits in memory.
Historical fits and the source of the current constants are retained under the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).

## Current contract

The worker floor is 2 GiB. The point, freeze-per-step and look-ahead native
paths each use a 2 MiB-per-spacecraft estimate; the offline surrogate adds no
native-GRAM per-spacecraft term. `SPACEAGORA_GRAM_SAT_MEMORY_MB` overrides the
native per-spacecraft estimate across paths. See `machine_topology.jl` for the
byte constants and admission predicates, and
`test/unit/parallel/gram_memory_model_tests.jl` for their relationships.

The native path table is explicit even while its entries are equal.
`OuterRouteFeatures.density_family` does not distinguish the environment-selected
point, frozen and look-ahead modes; the process-worker estimate currently uses
the point default at its call site. A future fit that separates those paths
would also need a way to carry the chosen path to that call site.

## Measurement method

`benchmarks/studies/gram_memory_footprint` launches a separate one-thread Julia
process for each density-path/constellation-size point. This matches the
process-pool worker shape and prevents an earlier point's heap from becoming a
later point's baseline. A one-spacecraft warm-up pays the first solver/native
allocations before the N-spacecraft configuration is measured.

The driver records peak RSS, live RSS before and after collection, and Linux
`VmRSS` when available. The recorded baseline RSS is not subtracted from peak
RSS; page release after collection is not synchronous. Fit an intercept and a
slope to the peak measurements, and include repeated solves because pool
workers serve more than one sample.

| Path | Configuration |
|---|---|
| `none` | no atmosphere control |
| `gram_point` | direct native point density, no freeze or vacuum cache |
| `gram_freeze_step` | density frozen per accepted step |
| `gram_lookahead` | vacuum-predicted look-ahead cache; the study extends its horizon and deviation limits past the mission |
| `gram_surrogate` | offline surrogate table |

The study uses the `paper_scenarios` geometry, force stack, tolerances and save
settings, and leaves state isolation at the production default. Record any
changed setting with the new fit.

```bash
julia --project=. benchmarks/studies/gram_memory_footprint/run_footprint.jl
julia --project=. benchmarks/studies/gram_memory_footprint/fit_footprint.jl benchmarks/studies/gram_memory_footprint/results/gram_memory_footprint_<host>.csv
julia --project=. --threads=8 benchmarks/studies/gram_memory_footprint/route_fit_table.jl
```

Study controls include `GMF_PATHS`, `GMF_SIZES`, `GMF_MISSION_S`, `GMF_GRAVITY`,
`GMF_ISOLATE`, `GMF_SOLVE_REPEATS`, `GMF_TIMEOUT_S`, `GMF_MIN_FREE_GB` and
`GMF_SUFFIX`. Keep the raw rows and host/revision metadata with the fit in the
private results archive rather than adding them to this source page.

## Limits

The fit does not price a density-service pool or isolated native pools within
one process, and it does not model inner-thread width. Long missions, denser
output and other force stacks need their own measurement. Surrogate queries
must use a table built for the scenario epoch; a surrogate is not a native
worker-memory observation. Timings gathered on a shared host are diagnostic,
not portable performance claims.
