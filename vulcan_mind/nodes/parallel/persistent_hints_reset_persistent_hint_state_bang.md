---
id: parallel.persistent_hints_reset_persistent_hint_state_bang
label: reset_persistent_hint_state!
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: reset_persistent_hint_state!
  lines:
  - 160
  - 160
inputs:
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
  description: Return value of `reset_persistent_hint_state!`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# reset_persistent_hint_state!

## Purpose
Discards the in-memory hint state so the next operation reloads from disk, used by tests and by profile switches that must not inherit learned allotments.

## Design & Implementation
Under `_persistent_hint_lock`, it clears `loaded` and `dirty`, blanks the path and empties the history dictionary. Clearing `loaded` is what makes the next `_ensure_persistent_hint_state_loaded!` re-read the file. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `reset_persistent_hint_state!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:171-171`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It does not touch the file, and because `dirty` is cleared any unsaved observations are lost rather than flushed; the `atexit` registration flag is not reset, which is correct for one process but means the hook keeps pointing at whichever path is loaded at exit.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 160.
