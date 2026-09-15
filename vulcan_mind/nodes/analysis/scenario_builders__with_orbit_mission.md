---
id: analysis.scenario_builders__with_orbit_mission
label: _with_orbit_mission
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _with_orbit_mission
  lines:
  - 453
  - 453
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: target_orbits
  type: Int
  units: n/a
  required: true
  description: Positional argument `target_orbits`.
- id: mission_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_time_s`.
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
  description: Return value of `_with_orbit_mission`.
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

# _with_orbit_mission

## Purpose
Switches a configuration's mission definition to orbit-counted termination with the given orbit count and a time bound.

## Design & Implementation
Rebuilds `MissionConfiguration` with `mission_type=MissionOrbits`, `number_of_orbits=target_orbits` and `mission_time=mission_time_s`, keeping `keplerian`, `orientation_sim`, `num_steps_to_save` and `data_rate` from the existing one, then rebuilds the `SimulationConfiguration` around it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `target_orbits` | Int | n/a | yes | Positional argument `target_orbits`. |
| in | `mission_time_s` | Float64 | n/a | yes | Positional argument `mission_time_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_with_orbit_mission`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:584-584`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:459-459`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:468-468`
<!-- vulcan:connections:end -->

## Limitations
The time bound is computed from the two-body period by the caller, so if drag shortens orbits substantially the run terminates on orbit count well before the time bound, and vice versa if the estimate is short.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 453.
