---
id: parallel.env_config_persistent_hints_exploration
label: persistent_hints_exploration
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: persistent_hints_exploration
  lines:
  - 83
  - 83
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
  type: Float64
  units: n/a
  description: Return value of `persistent_hints_exploration`.
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

# persistent_hints_exploration

## Purpose
Provides the exploration constant used by the persisted-hint selection rule (an upper-confidence-bound style bonus) so that thread counts with few samples are still occasionally tried.

## Design & Implementation
Parses `SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION` via `parse_float_env` with default `1.5`. If the parsed value `c` is not strictly positive the function returns `1.5` rather than throwing, so zero or negative exploration silently reverts to the default. Returns `Float64`; malformed text still raises `ArgumentError` from the parser.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `persistent_hints_exploration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:391-391`
- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:276-276`

**Downstream**

- `callees` → [[parallel.env_config_parse_float_env|parse_float_env]] · `callers` · call · `src/parallel/policy/env_config.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
`NaN` compares false in `c > 0.0` and therefore also falls back to 1.5 without any diagnostic. `Inf` passes the check and would make exploration dominate the score. The fallback-instead-of-error behaviour is inconsistent with `adaptive_delta` and `adaptive_rho`, which throw on out-of-range input.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 83.
