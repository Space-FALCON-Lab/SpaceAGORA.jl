---
id: vehicle.robotics_clotharmpose
label: ClothArmPose
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmPose
  lines:
  - 53
  - 53
inputs:
- id: base_position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `base_position`.
- id: base_quaternion
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `base_quaternion`.
- id: joint_origins
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `joint_origins`.
- id: joint_axes_world
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `joint_axes_world`.
- id: link_com_positions
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_com_positions`.
- id: link_tip_positions
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_tip_positions`.
- id: link_quaternions
  type: Vector{SVector{4, Float64}}
  units: n/a
  required: true
  description: Field `link_quaternions`.
- id: end_effector_position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `end_effector_position`.
- id: end_effector_quaternion
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `end_effector_quaternion`.
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
  type: ClothArmPose
  units: n/a
  description: Constructed `ClothArmPose`.
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

# ClothArmPose

## Purpose
Everything forward kinematics produces for one joint configuration: every joint origin and axis in world coordinates, each link's centre of mass, tip and orientation, and the end-effector frame.

## Design & Implementation
Immutable, echoing the base position and quaternion, then five parallel vectors indexed by link — `joint_origins`, `joint_axes_world`, `link_com_positions`, `link_tip_positions`, `link_quaternions` — and finally `end_effector_position` and `end_effector_quaternion`, which equal the last link's tip and orientation. Storing world-frame axes alongside origins is what lets a geometric Jacobian be assembled without re-running FK.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `base_position` | SVector{3, Float64} | n/a | yes | Field `base_position`. |
| in | `base_quaternion` | SVector{4, Float64} | n/a | yes | Field `base_quaternion`. |
| in | `joint_origins` | Vector{SVector{3, Float64}} | n/a | yes | Field `joint_origins`. |
| in | `joint_axes_world` | Vector{SVector{3, Float64}} | n/a | yes | Field `joint_axes_world`. |
| in | `link_com_positions` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_com_positions`. |
| in | `link_tip_positions` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_tip_positions`. |
| in | `link_quaternions` | Vector{SVector{4, Float64}} | n/a | yes | Field `link_quaternions`. |
| in | `end_effector_position` | SVector{3, Float64} | n/a | yes | Field `end_effector_position`. |
| in | `end_effector_quaternion` | SVector{4, Float64} | n/a | yes | Field `end_effector_quaternion`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmPose | n/a | — | Constructed `ClothArmPose`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:198-198`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The five per-link vectors are heap-allocated on every FK call, so the finite-difference Jacobian in `_ee_position_jacobian`, which calls FK `2n` times, allocates proportionally.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 53.
