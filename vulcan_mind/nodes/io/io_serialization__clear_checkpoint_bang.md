---
id: io.io_serialization__clear_checkpoint_bang
label: _clear_checkpoint!
kind: function
source:
  file: src/io/serialization/io_serialization.jl
  symbol: _clear_checkpoint!
  lines:
  - 83
  - 83
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
  type: Nothing
  units: n/a
  description: Return value of `_clear_checkpoint!`; mutates `args` in place. Returns
    `nothing`.
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

# _clear_checkpoint!

## Purpose
Removes both files of a checkpoint pair, used when a run completes successfully and its resume state should not be picked up by the next invocation.

## Design & Implementation
Resolves the data and manifest paths from `IOConfig._checkpoint_paths` and removes each with `rm` guarded by an `isfile` test and called with `force`, so neither a missing file nor a read-only flag raises. Deleting both is what keeps the pair consistent: a surviving manifest pointing at a deleted payload would otherwise look like a valid checkpoint to a tool that only inspects the manifest.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_clear_checkpoint!`; mutates `args` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/serialization/io_serialization.jl`

**Downstream**

- `callees` → [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:84-84`
- `callees` → [[simulation.resume_checkpoint__checkpoint_paths|_checkpoint_paths]] · `callers` · call · `src/io/serialization/io_serialization.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
The two removals are independent, so a failure between them leaves a half-deleted pair; nothing verifies afterwards that both are gone.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl` line 83.
