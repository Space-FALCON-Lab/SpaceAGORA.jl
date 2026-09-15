---
id: io.io_outputs__write_results_csv_bang
label: _write_results_csv!
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _write_results_csv!
  lines:
  - 96
  - 96
inputs:
- id: results_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `results_df`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: String
  units: n/a
  description: Return value of `_write_results_csv!`; mutates `results_df` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# _write_results_csv!

## Purpose
`_write_results_csv!` writes the results `DataFrame` to the run's primary CSV path atomically, and resolves the case where another process wrote the same file concurrently by preserving the other writer's output under a collision path before overwriting. It returns the path actually written.

## Design & Implementation
Signature `(results_df::DataFrame, args)::String`. It resolves `primary_path = IOConfig._results_csv_path(args)`, records `started_s = time()` and whether the file `existed_before`, then attempts `IOSerialization._atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=false)`. That call throws `ArgumentError` when the destination already exists. In the `catch`, if the error is an `ArgumentError` and the file is present, it reads `stat(primary_path).mtime` (falling back to `0.0` on failure) and treats the situation as a concurrent write when the file did not exist at start or its mtime is at or after `started_s`; in that case the frame is first written to `IOConfig._collision_results_csv_path(args)`, and then the primary is written with `force=true`. Any other error is rethrown.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results_df` | DataFrame | n/a | yes | Positional argument `results_df`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_write_results_csv!`; mutates `results_df` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`
- [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:127-127`
- [[simulation.persistence__find_sample_value|_find_sample_value]] · `callees` → `callers` · call · `src/simulation/engine/persistence.jl:21-21`

**Downstream**

- `callees` → [[io.io_config__collision_results_csv_path|_collision_results_csv_path]] · `callers` · call · `src/io/outputs/io_outputs.jl:111-111`
- `callees` → [[io.io_config__results_csv_path|_results_csv_path]] · `callers` · call · `src/io/outputs/io_outputs.jl:97-97`
- `callees` → [[io.io_serialization__atomic_write_file|_atomic_write_file]] · `callers` · call · `src/io/outputs/io_outputs.jl:101-101`
- `callees` → [[simulation.persistence__collision_results_csv_path|_collision_results_csv_path]] · `callers` · call · `src/io/outputs/io_outputs.jl:111-111`
- `callees` → [[simulation.persistence__results_csv_path|_results_csv_path]] · `callers` · call · `src/io/outputs/io_outputs.jl:97-97`
<!-- vulcan:connections:end -->

## Limitations
The concurrency heuristic relies on filesystem mtime resolution and clock agreement between `time()` and the file system; on coarse-mtime filesystems a pre-existing stale file can be misclassified and a fresh collision copy written. When a collision is detected the caller's own data, not the competing writer's, is written to the collision path, and the primary is then overwritten, so the other process's output is lost. The collision write itself uses `force=false` and can throw if that path also exists. The whole `DataFrame` is serialised twice on a collision.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 96.
