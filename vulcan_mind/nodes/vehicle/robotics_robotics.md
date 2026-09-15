---
id: vehicle.robotics_robotics
label: Robotics
kind: module
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: Robotics
  lines:
  - 2
  - 2
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
  description: Value produced by this symbol.
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

# Robotics

## Purpose
The kinematics module for the cloth-handling serial arm: model types, forward kinematics, damped least-squares inverse kinematics and surface target selection.

## Design & Implementation
Depends only on LinearAlgebra and StaticArrays. It defines six immutable types — base pose, link, joint, model, pose and state — and exports them alongside `default_cloth_arm_model`, the FK entry points, `cloth_ik`, `cloth_total_reach` and `closest_surface_target`. Quaternion arithmetic is kept private in `_unit_quat`, `_quat_mul`, `_quat_from_axis_angle` and `_rot`, all `@inline` on static vectors so a forward-kinematics pass allocates only the per-link output vectors. The scalar-last quaternion convention matches the rest of the vehicle model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Everything here is pure kinematics; joint velocities are carried in `ClothArmState` but never propagated to link velocities, which the dynamics layer computes separately. There is no collision checking between links.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 2.
