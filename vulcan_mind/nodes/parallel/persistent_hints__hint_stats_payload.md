---
id: parallel.persistent_hints__hint_stats_payload
label: _hint_stats_payload
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_stats_payload
  lines:
  - 9
  - 9
inputs:
- id: stats
  type: AdaptiveChoiceStats
  units: n/a
  required: true
  description: Positional argument `stats`.
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
  type: Dict{String,
  units: n/a
  description: Return value of `_hint_stats_payload`.
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

# _hint_stats_payload

## Purpose
Serialises one `AdaptiveChoiceStats` record into the string-keyed dictionary that the TOML hint file stores.

## Design & Implementation
Builds a `Dict{String,Any}` with the five fields — `samples`, `successes`, `failures`, `elapsed_sum_ns` and `elapsed_sq_sum_ns` — converting integers through `Int` and the elapsed sums through `Float64`, and clamping every value at zero. The clamps mean a corrupted in-memory record cannot write a negative count to disk that the loader would then reject.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `stats` | AdaptiveChoiceStats | n/a | yes | Positional argument `stats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{String, | n/a | — | Return value of `_hint_stats_payload`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints__save_persistent_hint_state_locked_bang|_save_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:126-126`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:14-14`
<!-- vulcan:connections:end -->

## Limitations
The elapsed sums are written as floating-point nanoseconds, so after very long histories the squared sum loses integer precision and the derived variance becomes approximate.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 9.
