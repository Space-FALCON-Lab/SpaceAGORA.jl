---
id: simulation.registry__simulation_engine_module
label: _simulation_engine_module
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: _simulation_engine_module
  lines:
  - 31
  - 31
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
  type: Any
  units: n/a
  description: Return value of `_simulation_engine_module`. Returns `getproperty(root,
    :SimulationEngine)`.
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

# _simulation_engine_module

## Purpose
Late-binding accessor that returns the `SimulationEngine` module from the package root, letting callback code (which is compiled before the engine module exists) reach engine functions for callback-stage sampling without a circular `using`.

## Design & Implementation
An `@inline` function with no arguments. It walks up from the module-level constant `_simulation_model_module` (the parent of the callbacks module, i.e. `SimulationModel`) to its own parent, the package root, and checks `isdefined(root, :SimulationEngine)`. If the module has not yet been defined it calls `error("SimulationEngine module is not available for callback stage sampling.")`; otherwise it returns `getproperty(root, :SimulationEngine)`. Returning a `Module` means the caller dereferences names dynamically, so this path is not type-stable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_simulation_engine_module`. Returns `getproperty(root, :SimulationEngine)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:19-19`
- [[simulation.event_callbacks_condition|condition]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:213-213`
- [[simulation.event_callbacks_condition_bang|condition!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:8-8`
- [[simulation.event_callbacks_get_drag_state_callback|get_drag_state_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:141-141`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:98-98`
- [[simulation.event_callbacks_get_orbit_end_callback|get_orbit_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:43-43`
- [[simulation.planet_frame_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:79-79`
- [[simulation.runtime__buffered_stage_environment_state|_buffered_stage_environment_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:74-74`
- [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:57-57`
- [[simulation.runtime__stage_environment_state|_stage_environment_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:207-207`
- [[simulation.save_fields__save_altitude|_save_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:91-91`
- [[simulation.save_fields__save_heat_load|_save_heat_load]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:146-146`
- [[simulation.save_fields__save_latitude_deg|_save_latitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:105-105`
- [[simulation.save_fields__save_longitude_deg|_save_longitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:119-119`
- [[simulation.save_fields__save_mass|_save_mass]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:155-155`
- [[simulation.save_fields__save_periapsis_altitude|_save_periapsis_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:77-77`
- [[simulation.save_fields__save_positions|_save_positions]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:20-20`
- [[simulation.save_fields__save_quaternion|_save_quaternion]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:163-163`
- [[simulation.save_fields__save_velocities|_save_velocities]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:28-28`
- [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:22-22`
- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:150-150`
- [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:8-8`
- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:79-79`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Relies on the exact module nesting `root -> SimulationModel -> Callbacks`; moving the callbacks module changes what `parentmodule` returns and breaks the lookup. `getproperty` on a `Module` is a dynamic lookup, so downstream calls through it incur runtime dispatch. It throws rather than degrading gracefully if the engine is not loaded.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 31.
