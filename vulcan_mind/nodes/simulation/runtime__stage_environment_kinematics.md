---
id: simulation.runtime__stage_environment_kinematics
label: _stage_environment_kinematics
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _stage_environment_kinematics
  lines:
  - 56
  - 56
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Any
  units: n/a
  description: Return value of `_stage_environment_kinematics`. Returns `(`.
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

# _stage_environment_kinematics

## Purpose
Gathers the geometric quantities the atmosphere models need — planet-fixed state, altitude, latitude and longitude — from one satellite's inertial state.

## Design & Implementation
Extracts position and velocity with `_extract_pos_vel`, then calls the engine's `sample_planet_frame` for satellite one at time `t` to obtain the inertial-to-planet rotation `l_pi`, planet-fixed position and velocity, and geodetic altitude, latitude and longitude. Returns them as a named tuple. Looking the engine module up through `_simulation_engine_module()` avoids a circular import between callbacks and engine.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_stage_environment_kinematics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__buffered_stage_environment_state|_buffered_stage_environment_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:73-73`
- [[simulation.runtime__stage_environment_state|_stage_environment_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:205-205`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:293-293`
- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:57-57`
- `callees` → [[simulation.runtime__extract_pos_vel|_extract_pos_vel]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:58-58`
<!-- vulcan:connections:end -->

## Limitations
`sample_planet_frame` is always called with satellite index one, so any per-satellite frame customisation is ignored here; the tuple copies every field, which is cheap for static vectors but adds up across many satellites per step.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 56.
