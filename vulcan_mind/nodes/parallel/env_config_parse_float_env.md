---
id: parallel.env_config_parse_float_env
label: parse_float_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: parse_float_env
  lines:
  - 33
  - 33
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Float64
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
  type: Float64
  units: n/a
  description: Return value of `parse_float_env`.
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

# parse_float_env

## Purpose
Reads a floating-point tuning parameter from the environment, returning the caller's `Float64` default when the variable is unset and raising a descriptive error on malformed text.

## Design & Implementation
`parse_float_env(name::String, default::Float64)::Float64`. The default is converted with `string(default)` so `get(ENV, name, ...)` always yields text, then `strip` removes whitespace. `parse(Float64, raw)` is wrapped in `try`/`catch`, and failures become `ArgumentError("<name> must be a float, got '<raw>'")`. No range checking is performed here; callers such as `adaptive_delta` and `adaptive_rho` enforce their own bounds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Float64 | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `parse_float_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.env_config_adaptive_delta|adaptive_delta]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:216-216`
- [[parallel.env_config_adaptive_rho|adaptive_rho]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:224-224`
- [[parallel.env_config_persistent_hints_exploration|persistent_hints_exploration]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:84-84`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:216-216`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Accepts `NaN`, `Inf`, and negative values unchanged, so every caller must validate. `string(default)` for values like `1.5` round-trips exactly but the pattern relies on Julia's shortest-representation printing. Empty environment strings throw rather than defaulting.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 33.
