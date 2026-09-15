---
id: io.io_config__collision_results_csv_path
label: _collision_results_csv_path
kind: function
source:
  file: src/io/config/io_config.jl
  symbol: _collision_results_csv_path
  lines:
  - 14
  - 14
inputs:
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
  description: Return value of `_collision_results_csv_path`.
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

# _collision_results_csv_path

## Purpose

Produces a results CSV path that cannot collide with another concurrently running simulation, for use when several processes or workers write into the same results directory.

## Design & Implementation

Formats the current UTC instant with `dateformat"yyyymmddTHHMMSSsss"` (millisecond resolution), then concatenates that stamp with `getpid()` and a `rand(UInt)` draw into a token, yielding `simulation_results.<stamp>.<pid>.<rand>.csv` under `args.simulation_settings.results_directory`. The three-part token defends against the three realistic clash sources: same millisecond, same host, same process.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_collision_results_csv_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:111-111`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

`rand(UInt)` draws from the task-local default RNG, so a run that reseeds the global RNG for reproducibility will make these tokens repeat across otherwise identical runs; uniqueness then rests on the timestamp and PID alone. PIDs are recycled by the operating system, and the returned path is never reserved, so a caller must still open the file with create semantics. The name is also not sortable back to a run identity beyond its timestamp.

## Provenance
Mapped from `src/io/config/io_config.jl` line 14.
