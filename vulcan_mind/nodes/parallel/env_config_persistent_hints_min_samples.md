---
id: parallel.env_config_persistent_hints_min_samples
label: persistent_hints_min_samples
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: persistent_hints_min_samples
  lines:
  - 88
  - 88
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
  type: Int
  units: n/a
  description: Return value of `persistent_hints_min_samples`.
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

# persistent_hints_min_samples

## Purpose
Sets how many timing observations a candidate thread count must accumulate before its persisted hint is trusted over exploration.

## Design & Implementation
Delegates to `parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES", 2)`, inheriting the clamp to at least 1 and the `ArgumentError` on non-integer text. Returns `Int`. The default of 2 means a single noisy measurement cannot lock in a choice.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `persistent_hints_min_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:392-392`
- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:277-277`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/parallel/policy/env_config.jl:89-89`
<!-- vulcan:connections:end -->

## Limitations
No upper bound, so a very large value effectively disables exploitation of persisted hints for short runs. The threshold is global across all telemetry buckets (`:density`, `:control`, `:multibody`, `:other`) with no per-bucket override.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 88.
