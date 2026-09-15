---
id: vehicle.geometry_properties_get_sc_area
label: get_SC_area
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_SC_area
  lines:
  - 87
  - 87
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: body
  type: Link
  units: n/a
  required: true
  description: Positional argument `body`.
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
  description: Return value of `get_SC_area`. Returns `get_SC_area(bodies)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# get_SC_area

## Purpose
Totals the reference area of the spacecraft bus in an assembly, the complement of the solar array total.

## Design & Implementation
Mirrors the solar array accessor exactly: a model method that traverses from `body` and a vector method that sums `ref_area` over links whose `root` flag is true. Root links are the box-shaped bus bodies in this convention. Having the two accessors as mirror images means the bus and array areas partition the assembly's links with no double counting and no gaps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_SC_area`. Returns `get_SC_area(bodies)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:197-197`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- [[vehicle.geometry_properties_get_sa_area|get_SA_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:83-83`

**Downstream**

- `callees` → [[vehicle.geometry_properties_get_normal_vector|get_normal_vector]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:108-108`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
A multi-root assembly sums every root body into one bus area, so a configuration with several independent root bodies cannot distinguish their contributions through this accessor.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl` line 87.
