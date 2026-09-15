---
id: vehicle.geometry_properties_get_tangent_vector
label: get_tangent_vector
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_tangent_vector
  lines:
  - 127
  - 127
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
  description: 'Return value of `get_tangent_vector`. Returns `normalized ? normalize(tangent)
    : tangent`.'
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

# get_tangent_vector

## Purpose
Expresses a body's own z-axis in the inertial frame, giving the in-plane direction that pairs with the normal for surface-force decomposition.

## Design & Implementation
Identical in structure to the normal accessor: fetch the body-to-inertial rotation from `rotate_to_inertial`, apply it to the static unit vector along z, and optionally `normalize` the result when the keyword is set. Choosing z as the tangent means the normal and tangent are orthogonal by construction in the body frame, so any deviation from orthogonality in the inertial result is attributable entirely to the rotation matrix.

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
| out | `result` | Any | n/a | — | Return value of `get_tangent_vector`. Returns `normalized ? normalize(tangent) : tangent`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- [[vehicle.geometry_properties_get_normal_vector|get_normal_vector]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:123-123`

**Downstream**

- `callees` → [[vehx.kinematics_rotate_to_inertial|rotate_to_inertial]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:129-129`
<!-- vulcan:connections:end -->

## Limitations
Only one tangent direction is returned, so a model needing a full surface basis must derive the third axis itself; the result is not re-orthogonalised against the normal, so a non-orthonormal rotation matrix propagates skew into both.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl` line 127.
