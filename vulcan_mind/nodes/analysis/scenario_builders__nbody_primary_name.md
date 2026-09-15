---
id: analysis.scenario_builders__nbody_primary_name
label: _nbody_primary_name
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _nbody_primary_name
  lines:
  - 42
  - 42
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  type: String
  units: n/a
  description: Return value of `_nbody_primary_name`.
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

# _nbody_primary_name

## Purpose
Converts a scenario planet name into the capitalised NAIF body name the N-body gravity model expects as its primary.

## Design & Implementation
Lowercases and strips the input and maps `earth`, `mars`, `venus`, `moon` and `titan` to their capitalised forms, raising `ArgumentError` otherwise. `@inline` with a `::String` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_nbody_primary_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:151-151`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The mapping is a literal branch, so it and `_planet_from_name` must be kept in step by hand; they already disagree on Titan.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 42.
