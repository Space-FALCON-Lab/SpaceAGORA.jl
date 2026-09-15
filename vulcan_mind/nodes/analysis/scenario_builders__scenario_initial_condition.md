---
id: analysis.scenario_builders__scenario_initial_condition
label: _scenario_initial_condition
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _scenario_initial_condition
  lines:
  - 525
  - 525
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: SimulationModel.AbstractInitialCondition
  units: n/a
  description: Return value of `_scenario_initial_condition`.
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

# _scenario_initial_condition

## Purpose
Produces the initial state for an orbit-events scenario, preferring an exact kernel-derived Cartesian state over published osculating elements.

## Design & Implementation
If the manifest supplies `initial_state_j2000_m`, it prints a notice and returns a `CartesianInitialCondition` from its six elements. Otherwise it forms periapsis radius from `Rp_e` plus `rp_altitude_m`, builds an `InitialCondition` from apoapsis radius, inclination, argument of periapsis, RAAN and true anomaly, and passes it through `_initial_condition_in_j2000` with the manifest's `element_frame`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.AbstractInitialCondition | n/a | — | Return value of `_scenario_initial_condition`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:557-557`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:546-546`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:531-531`
- `callees` → [[vehicle.model_cartesianinitialcondition|CartesianInitialCondition]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:532-532`
- `callees` → [[vehicle.model_initialcondition|InitialCondition]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:538-538`
<!-- vulcan:connections:end -->

## Limitations
When the Cartesian override is present the published elements are documentation only, and no consistency check between the two is performed.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 525.
