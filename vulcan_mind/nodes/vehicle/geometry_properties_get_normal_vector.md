---
id: vehicle.geometry_properties_get_normal_vector
label: get_normal_vector
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_normal_vector
  lines:
  - 112
  - 112
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
- id: root_index
  type: Int
  units: n/a
  required: true
  description: Positional argument `root_index`.
- id: normalized
  type: Any
  units: n/a
  required: false
  description: Keyword argument `normalized` (default `false`).
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
  description: 'Return value of `get_normal_vector`. Returns `normalized ? normalize(normal)
    : normal`.'
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

# get_normal_vector

## Purpose
Expresses a body's own x-axis in the inertial frame, which is the surface normal convention flat-plate aerodynamic and radiation models expect.

## Design & Implementation
Obtains the body-to-inertial rotation matrix from `rotate_to_inertial` for the given `body` and `root_index`, then applies it to the static unit vector along x. The `normalized` keyword defaults to false and, when set, passes the result through `normalize`, which matters because a rotation matrix assembled from a drifting quaternion is not exactly orthonormal and the product can be slightly off unit length.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `root_index` | Int | n/a | yes | Positional argument `root_index`. |
| in | `normalized` | Any | n/a | no | Keyword argument `normalized` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_normal_vector`. Returns `normalized ? normalize(normal) : normal`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- [[vehicle.geometry_properties_get_sc_area|get_SC_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:108-108`

**Downstream**

- `callees` → [[vehicle.geometry_properties_get_tangent_vector|get_tangent_vector]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:123-123`
- `callees` → [[vehx.kinematics_rotate_to_inertial|rotate_to_inertial]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:114-114`
<!-- vulcan:connections:end -->

## Limitations
The x-axis normal convention is implicit and shared with the `Facet` default rather than declared anywhere, so a body whose geometry was authored with a different normal convention is silently mis-oriented.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl` line 112.
