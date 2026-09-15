---
id: simulation.resume_checkpoint__clear_checkpoint_bang
label: _clear_checkpoint!
kind: function
source:
  file: src/simulation/engine/resume_checkpoint.jl
  symbol: _clear_checkpoint!
  lines:
  - 15
  - 15
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
  type: Any
  units: n/a
  description: Return value of `_clear_checkpoint!`; mutates `args` in place. Returns
    `SimulationModel.IOSerialization._clear_checkpoint!(args)`.
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

# _clear_checkpoint!

## Purpose
Deletes the resume checkpoint associated with a run so that a subsequent launch starts from the scenario's initial conditions instead of resuming from stale saved state.

## Design & Implementation
Implemented as `@inline _clear_checkpoint!(args) = SimulationModel.IOSerialization._clear_checkpoint!(args)`. The bang in the name reflects the filesystem side effect performed by the delegate, not mutation of `args` itself, which this frame leaves untouched. Keeping the alias here lets engine code call a short local name in the same namespace as `_load_checkpoint` and `_checkpoint_paths`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_clear_checkpoint!`; mutates `args` in place. Returns `SimulationModel.IOSerialization._clear_checkpoint!(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Deletion is irreversible and unguarded: there is no confirmation, no backup copy, and no return value inspected at this level to tell the caller whether a file was actually removed or was already absent. If another process holds the checkpoint file open, the platform-level failure propagates from IOSerialization, and partial cleanup of a multi-file checkpoint set is not rolled back.

## Provenance
Mapped from `src/simulation/engine/resume_checkpoint.jl` line 15.
