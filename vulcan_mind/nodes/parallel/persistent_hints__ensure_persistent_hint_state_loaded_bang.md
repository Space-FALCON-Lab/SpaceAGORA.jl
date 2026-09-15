---
id: parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang
label: _ensure_persistent_hint_state_loaded!
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _ensure_persistent_hint_state_loaded!
  lines:
  - 144
  - 144
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
  description: Return value of `_ensure_persistent_hint_state_loaded!`.
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

# _ensure_persistent_hint_state_loaded!

## Purpose
The entry point every hint operation calls first: load the history if needed and, once per process, register the save-at-exit hook.

## Design & Implementation
Takes `_persistent_hint_lock`, calls the locked loader, and if hints are enabled and `_persistent_hint_atexit_registered[]` is still false, sets the flag and registers an `atexit` closure that takes the same lock and saves. Registering inside the lock guarantees exactly one hook even under concurrent first calls.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_ensure_persistent_hint_state_loaded!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.persistent_hints__hint_record_observation_bang|_hint_record_observation!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:347-347`
- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:388-388`
- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:251-251`

**Downstream**

- `callees` → [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:146-146`
- `callees` → [[parallel.persistent_hints__save_persistent_hint_state_locked_bang|_save_persistent_hint_state_locked!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:151-151`
<!-- vulcan:connections:end -->

## Limitations
Saving only at exit means a process killed by a signal loses every observation since the last explicit save; there is no periodic flush.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 144.
