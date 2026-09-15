---
id: parallel.env_config_parse_nonnegative_int_env
label: parse_nonnegative_int_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: parse_nonnegative_int_env
  lines:
  - 43
  - 43
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Int
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `parse_nonnegative_int_env`.
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

# parse_nonnegative_int_env

## Purpose
Parses an integer environment knob that may legitimately be zero (for example a disabled trim budget) and clamps negative input up to 0.

## Design & Implementation
Mirrors `parse_thread_threshold_env` but ends with `max(0, value)` instead of `max(1, ...)`. Signature `parse_nonnegative_int_env(name::String, default::Int)::Int`; `parse(Int, strip(raw))` failures are rethrown as `ArgumentError("<name> must be an integer, got '<raw>'")`. Currently used by `adaptive_trim_quanta_budget` for `SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Int | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `parse_nonnegative_int_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.env_config__snapshot_auto_min_budget|_snapshot_auto_min_budget]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:210-210`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:210-210`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Negative values are silently coerced to 0 rather than rejected, so a typo such as `-1` disables the feature without feedback. Integer overflow of `parse(Int, ...)` on very long digit strings raises an `OverflowError` that is caught and reported as a generic integer error.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 43.
