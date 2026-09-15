---
id: simulation_a.assembly_get_callbacks
label: get_callbacks
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: get_callbacks
  lines:
  - 141
  - 194
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats_effectors_args
  type: Tuple{Int, Tuple, SimulationConfiguration}
  units: n/a
  required: true
  description: Spacecraft count, the dynamic effector tuple that decides which physics
    callbacks are required, and the run configuration.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: callback_set
  type: CallbackSet
  units: n/a
  description: Ordered `CallbackSet` handed to the ODE solver, containing every callback
    the configured mission actually needs.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# get_callbacks

## Purpose
`get_callbacks` is the single assembly point that turns a run configuration into the `CallbackSet` the integrator executes. It is the exported entry the simulation engine calls once per solve, and it decides which of the impact, planet-frame, density, thermal, orbit-end, entry-end, drag-state, navigation, control, guidance, quaternion-projection and data-saving callbacks are actually installed.

## Model & Assumptions
Installation is predicate-driven rather than unconditional: each `_requires_*` helper inspects the effector tuple and the mission configuration so that a vacuum orbit run does not pay for atmosphere or thermal callbacks. A second gate is `backbone_mode`, true when the engine solver policy is `:gravity_backbone_split`; in that mode only the impact callback and the orbit-end callback survive, because the split integration scheme evaluates the remaining physics on its own schedule. Save fields default to `default_save_fields(args)` when the caller passes none.

## Design & Implementation
Callbacks accumulate in an immutable tuple through the `_append_callback` and `_append_callbacks` helpers, which have `Nothing` and vector methods so a predicate that yields no callback is a no-op and effector families returning vectors splat cleanly. Keeping the accumulator a tuple preserves concrete element types, so `CallbackSet(callbacks...)` specialises rather than dispatching dynamically per step. `_callback_tolerances_for_phase`, defined alongside, supplies the per-phase absolute and relative tolerance templates used when a run switches between orbital and atmospheric error control.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats_effectors_args` | Tuple{Int, Tuple, SimulationConfiguration} | n/a | yes | Spacecraft count, the dynamic effector tuple that decides which physics callbacks are required, and the run configuration. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `callback_set` | CallbackSet | n/a | — | Ordered `CallbackSet` handed to the ODE solver, containing every callback the configured mission actually needs. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:210-210`

**Downstream**

- `callees` → [[simulation.assembly__append_callback|_append_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:161-161`
- `callees` → [[simulation.assembly__append_callbacks|_append_callbacks]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:181-181`
- `callees` → [[simulation.assembly__requires_drag_state_callback|_requires_drag_state_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:176-176`
- `callees` → [[simulation.assembly__requires_entry_end_callback|_requires_entry_end_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:172-172`
- `callees` → [[simulation.assembly__requires_orbit_end_callback|_requires_orbit_end_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:168-168`
- `callees` → [[simulation.assembly__requires_quaternion_projection_callback|_requires_quaternion_projection_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:185-185`
- `callees` → [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:160-160`
- `callees` → [[simulation.assembly__requires_thermal_callback|_requires_thermal_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:164-164`
- `callees` → [[simulation.event_callbacks_get_data_saving_callback|get_data_saving_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:189-189`
- `callees` → [[simulation.event_callbacks_get_drag_state_callback|get_drag_state_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:177-177`
- `callees` → [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:173-173`
- `callees` → [[simulation.event_callbacks_get_orbit_end_callback|get_orbit_end_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:169-169`
- `callees` → [[simulation.event_callbacks_get_quaternion_projection_callback|get_quaternion_projection_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:186-186`
- `callees` → [[simulation.navigation_guidance_callbacks_get_guidance_callbacks|get_guidance_callbacks]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:183-183`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:150-150`
- `callees` → [[simulation.save_fields__resolve_save_fields|_resolve_save_fields]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:149-149`
- `callees` → [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:182-182`
- `callees` → [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:152-152`
- `callees` → [[simulation_a.navigation_guidance_callbacks_get_navigation_callbacks|get_navigation_callbacks]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:181-181`
- `callees` → [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:156-156`
- `callees` → [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:161-161`
- `callees` → [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:165-165`
<!-- vulcan:connections:end -->

## Limitations
The resulting ordering is fixed by the source sequence; callbacks with side effects that depend on each other must be added in the correct place rather than relying on the solver. Duplicate save-field names raise `ArgumentError` through `_resolve_save_fields`. Extra callbacks supplied by the caller are appended last and are not validated against the built-in set.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl:140-194`.
