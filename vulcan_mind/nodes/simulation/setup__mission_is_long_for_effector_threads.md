---
id: simulation.setup__mission_is_long_for_effector_threads
label: _mission_is_long_for_effector_threads
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _mission_is_long_for_effector_threads
  lines:
  - 498
  - 498
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Bool
  units: n/a
  description: Return value of `_mission_is_long_for_effector_threads`.
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

# _mission_is_long_for_effector_threads

## Purpose
Classifies a mission as long or short for the purpose of deciding whether threading setup and cost-model warm-up will pay off over the run.

## Design & Implementation
Reads `args.mission_configuration`. When `mission_cfg.mission_type == SimulationModel.MissionOrbits` it returns `mission_cfg.number_of_orbits >= _effector_long_orbit_threshold()` (default 8); otherwise `mission_cfg.mission_time >= _effector_long_mission_threshold_s()` (default 5400 s). Both thresholds are re-read from the environment at call time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_mission_is_long_for_effector_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`

**Downstream**

- `callees` → [[simulation.setup__effector_long_mission_threshold_s|_effector_long_mission_threshold_s]] · `callers` · call · `src/simulation/engine/setup.jl:503-503`
- `callees` → [[simulation.setup__effector_long_orbit_threshold|_effector_long_orbit_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:501-501`
<!-- vulcan:connections:end -->

## Limitations
The two branches use incommensurable units, so an 8-orbit GEO mission and an 8-orbit LEO mission are both 'long' despite a 15× difference in wall time. Missions terminated by an event (altitude, impact) are judged by their nominal duration, not the actual one.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 498.
