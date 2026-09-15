---
id: parallel.persistent_hints__hint_workload_signature
label: _hint_workload_signature
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_workload_signature
  lines:
  - 186
  - 186
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: threshold
  type: Int
  units: n/a
  required: true
  description: Positional argument `threshold`.
- id: budget
  type: Int
  units: n/a
  required: true
  description: Positional argument `budget`.
- id: outer_active
  type: Bool
  units: n/a
  required: true
  description: Positional argument `outer_active`.
- id: heavy_only
  type: Bool
  units: n/a
  required: true
  description: Positional argument `heavy_only`.
- id: heavy_work
  type: Bool
  units: n/a
  required: true
  description: Positional argument `heavy_work`.
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
  type: String
  units: n/a
  description: Return value of `_hint_workload_signature`.
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

# _hint_workload_signature

## Purpose
Builds the string key under which threading observations are stored, identifying a workload shape well enough that history from one run applies to the next.

## Design & Implementation
Joins nine `key=value` tokens with `|`: the parallel profile and machine label from the `SPACEAGORA_PARALLEL_PROFILE` and `SPACEAGORA_PERF_MACHINE_LABEL` environment variables passed through `_safe_token`, the decision `source`, bucketed `num_items`, `threshold` and `budget`, and three boolean flags for outer parallelism, heavy-only and heavy-work. Including the machine label is what stops a laptop's history from steering a workstation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `threshold` | Int | n/a | yes | Positional argument `threshold`. |
| in | `budget` | Int | n/a | yes | Positional argument `budget`. |
| in | `outer_active` | Bool | n/a | yes | Positional argument `outer_active`. |
| in | `heavy_only` | Bool | n/a | yes | Positional argument `heavy_only`. |
| in | `heavy_work` | Bool | n/a | yes | Positional argument `heavy_work`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_hint_workload_signature`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:28-28`

**Downstream**

- `callees` → [[parallel.env_config__safe_token|_safe_token]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:195-195`
- `callees` → [[parallel.persistent_hints__hint_bucket|_hint_bucket]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:201-201`
<!-- vulcan:connections:end -->

## Limitations
The two environment variables are read on every call rather than once, and both default to `"default"`, so runs on different machines that forget to set the label pollute one shared history.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 186.
