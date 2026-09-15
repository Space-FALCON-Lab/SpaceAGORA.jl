---
id: parallel.types__persistenthintstate
label: _PersistentHintState
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: _PersistentHintState
  lines:
  - 101
  - 101
inputs:
- id: loaded
  type: Bool
  units: n/a
  required: false
  description: Field `loaded` (default `false`).
- id: dirty
  type: Bool
  units: n/a
  required: false
  description: Field `dirty` (default `false`).
- id: path
  type: String
  units: n/a
  required: false
  description: Field `path` (default `""`).
- id: history
  type: Dict{String, Dict{Int64, AdaptiveChoiceStats}}
  units: n/a
  required: false
  description: Field `history` (default `Dict{String, Dict{Int64, AdaptiveChoiceStats}}()`).
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
  type: _PersistentHintState
  units: n/a
  description: Constructed `_PersistentHintState` (keyword constructor via @kwdef).
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

# _PersistentHintState

## Purpose
Process-wide cache of persisted thread-allotment hints. It holds the on-disk path, the loaded history of `AdaptiveChoiceStats` keyed by workload signature and allotment, and the `loaded`/`dirty` flags that control lazy loading and at-exit flushing.

## Design & Implementation
`Base.@kwdef mutable struct` with `loaded = false`, `dirty = false`, `path = ""` and an empty `history::Dict{String, Dict{Int64, AdaptiveChoiceStats}}`. The single instance is `_persistent_hint_state::Ref{_PersistentHintState}` and every access is wrapped in `lock(_persistent_hint_lock)`; `adaptive_decision.jl` reads `loaded` and counts entries under that lock. `dirty` is set when a new observation updates `history` and cleared after a save; `_persistent_hint_atexit_registered::Ref{Bool}` ensures the flush hook is registered once.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `loaded` | Bool | n/a | no | Field `loaded` (default `false`). |
| in | `dirty` | Bool | n/a | no | Field `dirty` (default `false`). |
| in | `path` | String | n/a | no | Field `path` (default `""`). |
| in | `history` | Dict{String, Dict{Int64, AdaptiveChoiceStats}} | n/a | no | Field `history` (default `Dict{String, Dict{Int64, AdaptiveChoiceStats}}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentHintState | n/a | — | Constructed `_PersistentHintState` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The history is unbounded: signatures accumulate for every distinct `(source, num_items, threshold, budget, ...)` tuple seen and are never evicted. Because the outer `Dict{String, ...}` is keyed by a string signature, any change to `_hint_workload_signature` formatting orphans all persisted entries. The struct is not safe to read without the lock even though `Ref` access looks atomic, since `history` is a shared mutable dict.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 101.
