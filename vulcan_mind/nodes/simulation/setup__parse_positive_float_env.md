---
id: simulation.setup__parse_positive_float_env
label: _parse_positive_float_env
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _parse_positive_float_env
  lines:
  - 159
  - 159
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
  description: Return value of `_parse_positive_float_env`.
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

# _parse_positive_float_env

## Purpose
Parses a strictly positive `Float64` tuning parameter (time steps, nanosecond budgets, scale factors) from the engine environment, rejecting zero and negative values that would break divisions or sample counts.

## Design & Implementation
`_parse_positive_float_env(name::String, default::Float64)::Float64`. The default is stringified so `_engine_env_get(name, string(default))` always returns text, which is stripped and passed to `parse(Float64, raw)` inside `try`/`catch`; failures rethrow as `ArgumentError("<name> must be a floating-point value, got '<raw>'")`. A second guard `parsed > 0.0 || throw(ArgumentError("<name> must be > 0.0, got <parsed>"))` enforces positivity. Used by every `*_dt_s`, `*_ns_threshold`, and `*_scale` accessor in this file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Float64 | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_parse_positive_float_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:647-647`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__effector_cost_ns_per_item_default|_effector_cost_ns_per_item_default]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:418-418`
- [[simulation.setup__effector_long_mission_threshold_s|_effector_long_mission_threshold_s]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:470-470`
- [[simulation.setup__effector_outer_work_scale|_effector_outer_work_scale]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:466-466`
- [[simulation.setup__effector_work_ns_per_worker_threshold|_effector_work_ns_per_worker_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:462-462`
- [[simulation.setup__nbody_ephemeris_cache_dt_s|_nbody_ephemeris_cache_dt_s]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:209-209`
- [[simulation.setup__planet_frame_cache_dt_s|_planet_frame_cache_dt_s]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:221-221`
- [[simulation.setup__rhs_flat_cost_heterogeneity_threshold|_rhs_flat_cost_heterogeneity_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:793-793`
- [[simulation.setup__rhs_flat_packet_heterogeneity_threshold|_rhs_flat_packet_heterogeneity_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:446-446`
- [[simulation.setup__rhs_flat_packet_target_min_ns|_rhs_flat_packet_target_min_ns]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:430-430`
- [[simulation.setup__rhs_flat_packet_work_ns_threshold|_rhs_flat_packet_work_ns_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:442-442`
- [[simulation.setup__rhs_flat_work_ns_threshold|_rhs_flat_work_ns_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:785-785`
- [[simulation.setup__rhs_flat_work_per_worker_ns_threshold|_rhs_flat_work_per_worker_ns_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:789-789`
- [[simulation.setup__srp_ephemeris_cache_dt_s|_srp_ephemeris_cache_dt_s]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:197-197`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
`Inf` satisfies `> 0.0` and is accepted, which for `dt_s` accessors would yield a single-sample cache. `NaN` is rejected only because the comparison is false, with a message that prints `NaN`. Empty environment strings throw instead of falling back to the default.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 159.
