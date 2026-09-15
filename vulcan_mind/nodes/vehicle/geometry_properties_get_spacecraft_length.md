---
id: vehicle.geometry_properties_get_spacecraft_length
label: get_spacecraft_length
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_spacecraft_length
  lines:
  - 45
  - 45
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
  description: Return value of `get_spacecraft_length`. Returns `max_length`.
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

# get_spacecraft_length

## Purpose
Reports the longest characteristic dimension in an assembly, used as the reference length for aerodynamic and control scaling.

## Design & Implementation
Traverses from `body` through `traverse_bodies` to collect every link attached to it, then walks that collection tracking the largest value of `b.dims[1]`, starting from zero. Taking the first dimension of each link treats `dims` as ordered with the largest extent first, which is the convention the link constructors establish.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_spacecraft_length`. Returns `max_length`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/geometry_properties.jl`

**Downstream**

- `callees` → [[vehicle.geometry_properties_get_sa_area|get_SA_area]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:58-58`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:48-48`
<!-- vulcan:connections:end -->

## Limitations
Only `dims[1]` is examined, so a link whose longest extent is stored in another slot is under-reported; an assembly of links that all have zero or negative `dims[1]` returns the initial zero rather than signalling that no length was found.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl` line 45.
