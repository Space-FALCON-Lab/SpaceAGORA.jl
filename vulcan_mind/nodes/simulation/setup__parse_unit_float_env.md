---
id: simulation.setup__parse_unit_float_env
label: _parse_unit_float_env
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _parse_unit_float_env
  lines:
  - 170
  - 170
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
  description: Return value of `_parse_unit_float_env`.
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

# _parse_unit_float_env

## Purpose
Parses a fractional parameter constrained to the half-open interval (0, 1], used for exponential-moving-average weights and overhead ratios where 0 would freeze the estimator and values above 1 are meaningless.

## Design & Implementation
Same structure as `_parse_positive_float_env`: reads `_engine_env_get(name, string(default))`, strips, parses with `parse(Float64, ...)` under `try`/`catch` (rethrowing `ArgumentError("<name> must be a floating-point value, got '<raw>'")`), then requires `0.0 < parsed <= 1.0` or throws `ArgumentError("<name> must satisfy 0.0 < value <= 1.0, got <parsed>")`. Consumers are `_effector_cost_ema_alpha` (default 0.2) and `_rhs_flat_packet_overhead_disable_ratio` (default 0.10).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Float64 | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_parse_unit_float_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__effector_cost_ema_alpha|_effector_cost_ema_alpha]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:426-426`
- [[simulation.setup__rhs_flat_packet_overhead_disable_ratio|_rhs_flat_packet_overhead_disable_ratio]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:450-450`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:171-171`
<!-- vulcan:connections:end -->

## Limitations
Exactly 1.0 is allowed, which for an EMA alpha means no smoothing at all. The default is not validated against the same range at definition time, so a bad literal default would only fail when the variable is unset. Reads the environment on every call.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 170.
