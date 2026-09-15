---
id: analysis.scenario_builders__with_environment_wind
label: _with_environment_wind
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _with_environment_wind
  lines:
  - 367
  - 367
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: include_wind
  type: Bool
  units: n/a
  required: true
  description: Positional argument `include_wind`.
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
  type: SimulationConfiguration
  units: n/a
  description: Return value of `_with_environment_wind`.
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

# _with_environment_wind

## Purpose
Returns a copy of a configuration with the environment's wind flag set as the scenario requests.

## Design & Implementation
Rebuilds `EnvironmentModel` field by field from the existing one with `wind=include_wind`, then rebuilds `SimulationConfiguration` around it, passing every other section through unchanged. Both structs are immutable, so this copy-with-change pattern is the only way to alter one flag.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `include_wind` | Bool | n/a | yes | Positional argument `include_wind`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_with_environment_wind`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:585-585`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:624-624`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[core.simulation_configuration_environmentmodel|EnvironmentModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:369-369`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:380-380`
<!-- vulcan:connections:end -->

## Limitations
The explicit field lists mean a field added to either struct later is dropped here until the function is updated; `solver_config` in particular is not forwarded.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 367.
