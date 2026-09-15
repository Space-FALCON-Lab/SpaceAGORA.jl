---
id: parallel.policy_telemetry__record_policy_decision_bang
label: _record_policy_decision!
kind: function
source:
  file: src/parallel/policy/policy_telemetry.jl
  symbol: _record_policy_decision!
  lines:
  - 8
  - 8
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
- id: mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `mode`.
- id: threshold
  type: Int
  units: n/a
  required: true
  description: Positional argument `threshold`.
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
- id: adaptive_enabled
  type: Bool
  units: n/a
  required: true
  description: Positional argument `adaptive_enabled`.
- id: desire
  type: Int
  units: n/a
  required: true
  description: Positional argument `desire`.
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
- id: outer_active
  type: Bool
  units: n/a
  required: true
  description: Positional argument `outer_active`.
- id: allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Positional argument `allow_with_outer`.
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
- id: use_threads
  type: Bool
  units: n/a
  required: true
  description: Positional argument `use_threads`.
- id: signature
  type: String
  units: n/a
  required: true
  description: Positional argument `signature`.
- id: hint_allotment
  type: Int64
  units: n/a
  required: true
  description: Positional argument `hint_allotment`.
- id: hint_confidence
  type: Float64
  units: n/a
  required: true
  description: Positional argument `hint_confidence`.
- id: hint_regret_ns
  type: Float64
  units: n/a
  required: true
  description: Positional argument `hint_regret_ns`.
- id: hints_loaded
  type: Bool
  units: n/a
  required: true
  description: Positional argument `hints_loaded`.
- id: hints_entries
  type: Int64
  units: n/a
  required: true
  description: Positional argument `hints_entries`.
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
  description: Return value of `_record_policy_decision!`; mutates `source` in place.
    Returns `nothing`.
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

# _record_policy_decision!

## Purpose
Writes one threading decision into the active context's telemetry record, preserving both the running totals and the full detail of the most recent decision.

## Design & Implementation
Takes `_policy_telemetry_lock` for the whole body. It increments `decisions_total`, and conditionally `threads_enabled_total`, `policy_threading_proposed_total` and `adaptive_decisions_total`. `_telemetry_bucket(source)` classifies the source into `:density`, `:control`, `:multibody` or a catch-all, and the matching pair of per-bucket counters is bumped. It then stores every argument into the `last_*` fields, clamping the numeric ones with `max` so a threshold, budget, desire or allotment recorded as zero or negative reads back as at least one and confidences and regrets never read back negative.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `mode` | Symbol | n/a | yes | Positional argument `mode`. |
| in | `threshold` | Int | n/a | yes | Positional argument `threshold`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `budget` | Int | n/a | yes | Positional argument `budget`. |
| in | `adaptive_enabled` | Bool | n/a | yes | Positional argument `adaptive_enabled`. |
| in | `desire` | Int | n/a | yes | Positional argument `desire`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `outer_active` | Bool | n/a | yes | Positional argument `outer_active`. |
| in | `allow_with_outer` | Bool | n/a | yes | Positional argument `allow_with_outer`. |
| in | `heavy_only` | Bool | n/a | yes | Positional argument `heavy_only`. |
| in | `heavy_work` | Bool | n/a | yes | Positional argument `heavy_work`. |
| in | `use_threads` | Bool | n/a | yes | Positional argument `use_threads`. |
| in | `signature` | String | n/a | yes | Positional argument `signature`. |
| in | `hint_allotment` | Int64 | n/a | yes | Positional argument `hint_allotment`. |
| in | `hint_confidence` | Float64 | n/a | yes | Positional argument `hint_confidence`. |
| in | `hint_regret_ns` | Float64 | n/a | yes | Positional argument `hint_regret_ns`. |
| in | `hints_loaded` | Bool | n/a | yes | Positional argument `hints_loaded`. |
| in | `hints_entries` | Int64 | n/a | yes | Positional argument `hints_entries`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_record_policy_decision!`; mutates `source` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:104-104`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:30-30`
- `callees` → [[parallel.env_config__telemetry_bucket|_telemetry_bucket]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations
Only the single most recent decision is retained in the `last_*` fields, so a burst of decisions leaves no per-decision history; the `max` clamps quietly repair out-of-range inputs rather than reporting that the caller produced them.

## Provenance
Mapped from `src/parallel/policy/policy_telemetry.jl` line 8.
