---
id: simulation.from_env__parse_bool
label: _parse_bool
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_bool
  lines:
  - 5
  - 5
inputs:
- id: raw
  type: Any
  units: n/a
  required: true
  description: Positional argument `raw`.
- id: default
  type: Bool
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
  type: Any
  units: n/a
  description: Return value of `_parse_bool`. Returns `default`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _parse_bool

## Purpose
Lenient boolean parser for environment values that accepts the common truthy and falsy spellings and returns a caller default for anything else, including an unset variable.

## Design & Implementation
`_parse_bool(raw, default::Bool)`. If `raw === nothing` (the sentinel returned by `get(env, name, nothing)`) it returns `default`. Otherwise it normalises with `lowercase(strip(String(raw)))` and returns `true` for `"1"`, `"true"`, `"yes"`, `"on"`, `false` for `"0"`, `"false"`, `"no"`, `"off"`, and `default` for any other token. Used throughout `simulation_engine_config_from_env` for the parallel, runtime-policy and artifact flags.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | Any | n/a | yes | Positional argument `raw`. |
| in | `default` | Bool | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_parse_bool`. Returns `default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:164-164`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Malformed values are silently coerced to the default rather than raising, so a typo such as `SPACEAGORA_SAVE_BUNDLE=flase` is indistinguishable from leaving it unset. `String(raw)` on a non-string, non-`nothing` value (an `Int` 1, say) throws a `MethodError`. The accepted token set is duplicated from `ParallelPolicy.parse_bool_env` rather than shared.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 5.
