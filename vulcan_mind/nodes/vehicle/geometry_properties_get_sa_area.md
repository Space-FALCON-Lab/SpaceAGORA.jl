---
id: vehicle.geometry_properties_get_sa_area
label: get_SA_area
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_SA_area
  lines:
  - 62
  - 62
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
  description: Return value of `get_SA_area`. Returns `get_SA_area(bodies)`.
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

# get_SA_area

## Purpose
Totals the reference area of the solar array surfaces in an assembly, separating them from the bus for drag and radiation bookkeeping.

## Design & Implementation
Two methods. The model method traverses from `body` with `traverse_bodies` and delegates; the vector method accumulates `b.ref_area` over every link whose `root` flag is false. Non-root links are the flat plates in this model's convention, so the `root` flag is what distinguishes an array panel from the bus box without needing a separate component type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_SA_area`. Returns `get_SA_area(bodies)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:197-197`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- [[vehicle.geometry_properties_get_spacecraft_length|get_spacecraft_length]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:58-58`

**Downstream**

- `callees` → [[vehicle.geometry_properties_get_sc_area|get_SC_area]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:83-83`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
The non-root test is a proxy: any non-root link, whatever it physically represents, is counted as solar array area, so an appendage such as an antenna boom inflates this total.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl` line 62.
