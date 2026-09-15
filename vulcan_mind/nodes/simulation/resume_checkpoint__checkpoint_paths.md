---
id: simulation.resume_checkpoint__checkpoint_paths
label: _checkpoint_paths
kind: function
source:
  file: src/simulation/engine/resume_checkpoint.jl
  symbol: _checkpoint_paths
  lines:
  - 2
  - 2
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
  type: SimulationModel.IOConfig._checkpoint_paths
  units: n/a
  description: Return value of `_checkpoint_paths`. Returns `SimulationModel.IOConfig._checkpoint_paths(args)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _checkpoint_paths

## Purpose
Produces the concrete file paths used by the resume-checkpoint machinery for a run described by `args`, so writers and loaders agree on where state, metadata, and temporary files live.

## Design & Implementation
Defined as `@inline _checkpoint_paths(args) = SimulationModel.IOConfig._checkpoint_paths(args)`. Like its sibling `_checkpoint_directory`, it is a one-line re-export that keeps the naming convention owned by `SimulationModel.IOConfig`. Path construction, extension choice, and any run-identifier interpolation are all decided downstream; this frame contributes only the inlined call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOConfig._checkpoint_paths | n/a | — | Return value of `_checkpoint_paths`. Returns `SimulationModel.IOConfig._checkpoint_paths(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_serialization__clear_checkpoint_bang|_clear_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:84-84`
- [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:63-63`
- [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:35-35`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No existence, permission, or collision checking is performed, so two concurrently running simulations that resolve to the same run identifier will resolve to the same paths and can overwrite each other's checkpoints. Errors from malformed IO configuration propagate untranslated from IOConfig, and the returned paths are computed fresh on every call rather than cached.

## Provenance
Mapped from `src/simulation/engine/resume_checkpoint.jl` line 2.
