---
id: simulation.from_env__parse_float_opt
label: _parse_float_opt
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_float_opt
  lines:
  - 57
  - 57
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
- id: env_name
  type: String
  units: n/a
  required: true
  description: Positional argument `env_name`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_parse_float_opt`.
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

# _parse_float_opt

## Purpose
Parses an optional positive floating-point solver knob such as `SPACEAGORA_SYMPLECTIC_DT_S`, returning `nothing` for an empty string and raising a descriptive error for anything that is not a positive number.

## Design & Implementation
`_parse_float_opt(raw::String, env_name::String)::Union{Nothing, Float64}`. It strips whitespace, returns `nothing` when the result is empty, then `tryparse(Float64, s)`. A parse failure or a value `<= 0.0` throws `ArgumentError("<env_name> must be a positive number, got '<s>'.")`. Used for the symplectic, gravity-backbone and multirate slow time steps, each wrapped in `_parse_or_default` so strictness can be relaxed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `env_name` | String | n/a | yes | Positional argument `env_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_parse_float_opt`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:96-96`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Inf` parses successfully and passes the positivity check, and `NaN` fails the `> 0.0` test with the same generic message. Units (seconds) are implied by the variable name only. Zero is rejected even though some callers might want it to mean "auto".

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 57.
