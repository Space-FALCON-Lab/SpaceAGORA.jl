---
id: simulation_a.save_fields_default_save_fields
label: default_save_fields
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: default_save_fields
  lines:
  - 170
  - 191
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Run configuration supplying the spacecraft vector, whose length fixes
    the per-satellite column count, and the orientation-simulation flag.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: save_fields
  type: Vector{SaveField}
  units: n/a
  description: Ordered save-field descriptors, each pairing an output name with a
    getter closure and per-satellite column metadata, consumed by the data-saving
    callback.
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
# default_save_fields

## Purpose
`default_save_fields` defines what a SpaceAGORA run records when the caller does not specify an output schema. It is one of the three symbols `SimulationCallbacks` exports, and it establishes the default column layout that downstream analysis, plotting and regression comparisons rely on.

## Model & Assumptions
The default schema spans thirteen quantities: inertial position and velocity, altitude, geodetic latitude and longitude in degrees, mass, wind, the drag, lift and cross aerodynamic components, periapsis altitude, and instantaneous heat rate together with accumulated heat load. Every one is marked `per_satellite=true`, so the writer expands each into one column group per spacecraft. Attitude is conditional: a quaternion field is appended only when `args.mission_configuration.orientation_sim` is set, because the quaternion getter raises `ArgumentError` if the state carries no orientation component.

## Design & Implementation
Each entry is a `SaveField` pairing a symbolic name with a closure over the spacecraft count, so the count is resolved once at construction rather than re-derived on every solver step. Position, velocity and quaternion pass an explicit `column_prefix` of `pos`, `vel` and `q` so the emitted headers stay short and stable; the remaining fields take their name as the prefix. The companion `_resolve_save_fields` returns this default when the caller passes nothing and otherwise collects the caller's iterable, rejecting duplicate names with an `ArgumentError` in both cases. `_save_snapshot` then evaluates every getter into a `SaveData` dictionary at the saving callback's firing instants.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `args` | SimulationConfiguration | n/a | yes | Run configuration supplying the spacecraft vector, whose length fixes the per-satellite column count, and the orientation-simulation flag. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `save_fields` | Vector{SaveField} | n/a | — | Ordered save-field descriptors, each pairing an output name with a getter closure and per-satellite column metadata, consumed by the data-saving callback. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.save_fields__resolve_save_fields|_resolve_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:194-194`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:206-206`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:188-188`
- `callees` → [[simulation.save_fields__save_altitude|_save_altitude]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:175-175`
- `callees` → [[simulation.save_fields__save_cross|_save_cross]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:182-182`
- `callees` → [[simulation.save_fields__save_drag|_save_drag]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:180-180`
- `callees` → [[simulation.save_fields__save_heat_load|_save_heat_load]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:185-185`
- `callees` → [[simulation.save_fields__save_heat_rate|_save_heat_rate]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:184-184`
- `callees` → [[simulation.save_fields__save_latitude_deg|_save_latitude_deg]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:176-176`
- `callees` → [[simulation.save_fields__save_lift|_save_lift]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:181-181`
- `callees` → [[simulation.save_fields__save_longitude_deg|_save_longitude_deg]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:177-177`
- `callees` → [[simulation.save_fields__save_mass|_save_mass]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:178-178`
- `callees` → [[simulation.save_fields__save_periapsis_altitude|_save_periapsis_altitude]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:183-183`
- `callees` → [[simulation.save_fields__save_positions|_save_positions]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:173-173`
- `callees` → [[simulation.save_fields__save_quaternion|_save_quaternion]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:188-188`
- `callees` → [[simulation.save_fields__save_velocities|_save_velocities]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:174-174`
- `callees` → [[simulation.save_fields__save_wind|_save_wind]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:179-179`
- `callees` → [[simulation.save_fields_savefield|SaveField]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:173-173`
<!-- vulcan:connections:end -->

## Limitations
Getters read the shared buffers the density and thermal callbacks populate, so a field such as wind or heat rate is only meaningful when the callback that fills its buffer is installed for the run. Adding fields increases per-step allocation because most getters build a fresh vector per call. The schema is fixed at setup; fields cannot be added or removed once the solve has begun.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl:169-191`.
