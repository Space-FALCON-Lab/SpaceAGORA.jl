---
id: analysis.scenario_builders__initial_condition_in_j2000
label: _initial_condition_in_j2000
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _initial_condition_in_j2000
  lines:
  - 499
  - 499
inputs:
- id: ic
  type: InitialCondition
  units: n/a
  required: true
  description: Positional argument `ic`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: element_frame
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `element_frame`.
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
  description: Return value of `_initial_condition_in_j2000`.
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

# _initial_condition_in_j2000

## Purpose
Converts orbital elements published relative to a body's mean equator into the J2000 Cartesian state the propagator integrates in.

## Design & Implementation
Returns the elements untouched for `:j2000` and raises for any frame other than `:body_equator_inertial`. Otherwise it obtains ephemeris time from a `SpiceEphemeridesModel`, takes the planet-fixed rotation `l_pi` at that time, extracts the pole direction from its third row, builds the body-equator basis with `_body_equator_frame_rotation`, converts the elements to position and velocity with `orbitalelemtorv`, and rotates both into J2000 as a `CartesianInitialCondition`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ic` | InitialCondition | n/a | yes | Positional argument `ic`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `element_frame` | Symbol | n/a | yes | Positional argument `element_frame`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.AbstractInitialCondition | n/a | — | Return value of `_initial_condition_in_j2000`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:546-546`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:514-514`
- `callees` → [[envana.ana_scenario_builders_body_equator_frame_rotation|_body_equator_frame_rotation]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:513-513`
- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:510-510`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:511-511`
- `callees` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:509-509`
- `callees` → [[vehicle.model_cartesianinitialcondition|CartesianInitialCondition]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:515-515`
<!-- vulcan:connections:end -->

## Limitations
It uses the same pole model as propagation, which is the point, but that means a scenario whose published elements were referenced to a different pole epoch is silently misrotated by the pole drift between epochs.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 499.
