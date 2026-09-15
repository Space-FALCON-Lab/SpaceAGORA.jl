---
id: parallel.adaptive_decision_use_threads_policy
label: use_threads_policy
kind: function
source:
  file: src/parallel/policy/adaptive_decision.jl
  symbol: use_threads_policy
  lines:
  - 142
  - 142
inputs:
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: mode
  type: Symbol
  units: n/a
  required: true
  description: Keyword argument `mode`.
- id: threshold
  type: Int
  units: n/a
  required: true
  description: Keyword argument `threshold`.
- id: heavy_work
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `heavy_work` (default `true`).
- id: heavy_only
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `heavy_only` (default `false`).
- id: outer_active
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `outer_active` (default `false`).
- id: allow_with_outer
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `allow_with_outer` (default `true`).
- id: source
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `source` (default `:other`).
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
  type: Bool
  units: n/a
  description: Return value of `use_threads_policy`.
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

# use_threads_policy

## Purpose
Answers the single question a call site cares about — should this loop run threaded — by reducing the full threading policy decision to its boolean verdict.

## Design & Implementation
Forwards `num_items` and every keyword (`mode`, `threshold`, `heavy_work`, `heavy_only`, `outer_active`, `allow_with_outer`, `source`) unchanged to `thread_policy_decision` and returns only its `.use_threads` field, discarding the rest of the decision record. Declared `@inline` with a `::Bool` return annotation so it compiles away at the call site.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `mode` | Symbol | n/a | yes | Keyword argument `mode`. |
| in | `threshold` | Int | n/a | yes | Keyword argument `threshold`. |
| in | `heavy_work` | Bool | n/a | no | Keyword argument `heavy_work` (default `true`). |
| in | `heavy_only` | Bool | n/a | no | Keyword argument `heavy_only` (default `false`). |
| in | `outer_active` | Bool | n/a | no | Keyword argument `outer_active` (default `false`). |
| in | `allow_with_outer` | Bool | n/a | no | Keyword argument `allow_with_outer` (default `true`). |
| in | `source` | Symbol | n/a | no | Keyword argument `source` (default `:other`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `use_threads_policy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/adaptive_decision.jl`

**Downstream**

- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:152-152`
<!-- vulcan:connections:end -->

## Limitations
Discarding the decision record loses the reason the policy chose as it did, so a call site using this wrapper cannot report why threading was declined; use `thread_policy_decision` directly when that matters.

## Provenance
Mapped from `src/parallel/policy/adaptive_decision.jl` line 142.
