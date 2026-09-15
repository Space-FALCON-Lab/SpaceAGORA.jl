---
id: simulation.setup__effector_long_orbit_threshold
label: _effector_long_orbit_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_long_orbit_threshold
  lines:
  - 473
  - 473
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_effector_long_orbit_threshold`.
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

# _effector_long_orbit_threshold

## Purpose
Orbit count above which an orbit-based mission is considered long enough to justify effector threading overheads.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_LONG_ORBIT_THRESHOLD", 8)`, clamped to at least 1. `_mission_is_long_for_effector_threads` compares `mission_cfg.number_of_orbits >= threshold` when `mission_type == SimulationModel.MissionOrbits`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_effector_long_orbit_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__mission_is_long_for_effector_threads|_mission_is_long_for_effector_threads]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:501-501`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:474-474`
<!-- vulcan:connections:end -->

## Limitations
Eight orbits at a high altitude can be far longer in wall time than 5400 s, so the two long-mission criteria are not mutually consistent. The value is read from `ENV` each call rather than captured in `RhsPlanEnvConfig`.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 473.
