---
id: parallel.persistent_hints__hint_entry_count
label: _hint_entry_count
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_entry_count
  lines:
  - 1
  - 1
inputs:
- id: state
  type: _PersistentHintState
  units: n/a
  required: true
  description: Positional argument `state`.
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
  type: Int
  units: n/a
  description: Return value of `_hint_entry_count`.
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

# _hint_entry_count

## Purpose
Counts how many signature-and-allotment pairs the persistent hint history holds, reported in telemetry so users can see how much learned state a run started with.

## Design & Implementation
Iterates the values of `state.history`, each a dictionary keyed by allotment, and sums their lengths. Declared `@inline` with an `::Int` return. It is a pure read and takes no lock, relying on the caller to hold `_persistent_hint_lock`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | _PersistentHintState | n/a | yes | Positional argument `state`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_hint_entry_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:43-43`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:138-138`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It counts entries regardless of whether their statistics carry any samples, so entries created by `get!` during a lookup that never recorded an observation still count.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 1.
