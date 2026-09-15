---
id: parallel.persistent_hints__save_persistent_hint_state_locked_bang
label: _save_persistent_hint_state_locked!
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _save_persistent_hint_state_locked!
  lines:
  - 111
  - 111
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
  description: Return value of `_save_persistent_hint_state_locked!`.
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

# _save_persistent_hint_state_locked!

## Purpose
Writes the in-memory hint history back to disk atomically, but only when it was loaded, has changed, and persistence is enabled.

## Design & Implementation
Three guards short-circuit the common no-op cases. Rows are emitted in sorted signature and sorted allotment order, skipping empty buckets and zero-sample stats, so the file is deterministic and diffs cleanly across runs. The payload carries `schema_version` 1. It writes to `path * ".tmp"`, then `mv` with `force` over the real file, so a crash mid-write leaves the previous file intact. Finally `state.dirty` is cleared.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_save_persistent_hint_state_locked!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang|_ensure_persistent_hint_state_loaded!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:151-151`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:123-123`
- `callees` → [[parallel.env_config_persistent_hints_persist_enabled|persistent_hints_persist_enabled]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:113-113`
- `callees` → [[parallel.persistent_hints__hint_stats_payload|_hint_stats_payload]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:126-126`
<!-- vulcan:connections:end -->

## Limitations
The temporary filename is fixed, so two processes saving to the same hint path can clobber each other's temporary file; unlike the checkpoint writer in `IOSerialization`, no pid or random suffix is used.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 111.
