---
id: simulation.from_env__env_bool
label: _env_bool
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _env_bool
  lines:
  - 1
  - 1
inputs:
- id: v
  type: Bool
  units: n/a
  required: true
  description: Positional argument `v`.
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
  description: 'Return value of `_env_bool`. Returns `v ? "1" : "0"`.'
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

# _env_bool

## Purpose
Renders a Julia `Bool` as the `"1"`/`"0"` string convention used by every `SPACEAGORA_*` boolean environment variable, for use when a config is written back into an override dictionary.

## Design & Implementation
Defined as `@inline _env_bool(v::Bool) = v ? "1" : "0"`. It is the inverse direction of `_parse_bool`, which accepts `1/true/yes/on` and `0/false/no/off`; emitting only `"1"`/`"0"` keeps the round trip canonical. Used a dozen times inside `_engine_env_overrides`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Bool | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_env_bool`. Returns `v ? "1" : "0"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__engine_env_overrides|_engine_env_overrides]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:238-238`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only `Bool` is accepted; passing an integer or `nothing` is a `MethodError`. The chosen spelling is fixed, so any external tool that expects `true`/`false` must parse the numeric form.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 1.
