---
id: analysis.manifest_parsing__parse_positive_float_env
label: _parse_positive_float_env
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_positive_float_env
  lines:
  - 23
  - 23
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  description: Return value of `_parse_positive_float_env`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _parse_positive_float_env

## Purpose
Reads an optional positive float from an environment variable, distinguishing unset from present.

## Design & Implementation
Returns `nothing` when the variable is empty, otherwise parses as `Float64`, converting a parse failure into an `ArgumentError`, and rejects values at or below zero. `@inline` with a `Union{Nothing, Float64}` return. Used by the study tolerance override variables.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_parse_positive_float_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:647-647`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
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

- *none*
<!-- vulcan:connections:end -->

## Limitations
Zero is rejected even though a tolerance of zero could be meaningful for a strict comparison; callers cannot express that through the environment.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 23.
