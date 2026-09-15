---
id: io.io_serialization__atomic_write_file
label: _atomic_write_file
kind: function
source:
  file: src/io/serialization/io_serialization.jl
  symbol: _atomic_write_file
  lines:
  - 10
  - 10
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
- id: writer
  type: Function
  units: n/a
  required: true
  description: Positional argument `writer`.
- id: force
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `force` (default `true`).
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
  type: Any
  units: n/a
  description: Return value of `_atomic_write_file`. Returns `path`.
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

# _atomic_write_file

## Purpose
Writes a file such that a reader never observes a partial result, which is what makes a checkpoint safe to take while a run may be killed at any moment.

## Design & Implementation
Creates the destination directory, then builds a temporary name in that same directory from a leading dot, the target basename, the process id, the thread id and a random `UInt`, so concurrent writers on any thread of any process cannot collide. It invokes `writer` on the temporary path and moves it onto the target with `force`, the move being atomic because source and destination share a filesystem. A `finally` block removes the temporary if it still exists, so a throwing writer leaves no debris. Returns the destination path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `writer` | Function | n/a | yes | Positional argument `writer`. |
| in | `force` | Bool | n/a | no | Keyword argument `force` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_atomic_write_file`. Returns `path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:101-101`
- [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:131-131`
- [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:43-43`
- [[simulation.persistence__collision_results_csv_path|_collision_results_csv_path]] · `callees` → `callers` · call · `src/simulation/engine/persistence.jl:5-5`
- [[simulation.setup__write_nbody_ephemeris_cache_file_bang|_write_nbody_ephemeris_cache_file!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1620-1620`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Atomicity depends on the rename being within one filesystem; a destination directory that is a mount point or a symlink onto another device degrades the move to a copy and loses the guarantee. The directory entry is renamed but never fsynced, so a power loss can still lose the write.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl` line 10.
