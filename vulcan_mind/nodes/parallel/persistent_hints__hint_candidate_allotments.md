---
id: parallel.persistent_hints__hint_candidate_allotments
label: _hint_candidate_allotments
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_candidate_allotments
  lines:
  - 210
  - 210
inputs:
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: budget
  type: Int
  units: n/a
  required: true
  description: Positional argument `budget`.
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
  type: Vector{Int64}
  units: n/a
  description: Return value of `_hint_candidate_allotments`.
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

# _hint_candidate_allotments

## Purpose
Enumerates the thread allotments the hint bandit is allowed to choose between for a given item count and budget.

## Design & Implementation
Caps candidates at the smaller of `num_items` and `budget`, each at least one. It always includes 1, then 2, 4 and 8 where they fit under the cap, then half the cap rounded up and the cap itself, and returns the sorted, deduplicated `Int64` vector. Powers of two plus the cap and its half give a small, well-spread set to explore.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `budget` | Int | n/a | yes | Positional argument `budget`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Int64} | n/a | — | Return value of `_hint_candidate_allotments`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:37-37`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:214-214`
<!-- vulcan:connections:end -->

## Limitations
Between 8 and a large cap only the half and full values are tried, so on a 64-thread machine allotments of 16 or 24 are never evaluated even if they would be best.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 210.
