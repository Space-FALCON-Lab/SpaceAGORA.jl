---
id: vehicle.robotics__ee_position_jacobian
label: _ee_position_jacobian
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _ee_position_jacobian
  lines:
  - 238
  - 238
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: base_pose
  type: ClothArmBasePose
  units: n/a
  required: true
  description: Positional argument `base_pose`.
- id: q
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `q`.
- id: ee
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `ee`.
- id: h
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `h` (default `1.0e-6`).
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
  description: Return value of `_ee_position_jacobian`. Returns `J`.
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

# _ee_position_jacobian

## Purpose
Estimates how the end-effector position changes with each joint angle, the matrix the inverse-kinematics solver inverts.

## Design & Implementation
For each joint it perturbs the angle by a step of `max(h, h * |q[i]|)` with `h` defaulting to 1e-6, evaluates `cloth_fk` at plus and minus the step, and stores the central difference in column `i` of a three-by-`n` matrix, restoring the trial vector after each column. The `ee` argument is accepted but unused, since central differences need no reference evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q` | Vector{Float64} | n/a | yes | Positional argument `q`. |
| in | `ee` | SVector{3, Float64} | n/a | yes | Positional argument `ee`. |
| in | `h` | Float64 | n/a | no | Keyword argument `h` (default `1.0e-6`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_ee_position_jacobian`. Returns `J`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_cloth_ik|cloth_ik]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:277-277`

**Downstream**

- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/vehicle/robotics/robotics.jl:251-251`
<!-- vulcan:connections:end -->

## Limitations
Central differences cost two full FK passes per joint, each allocating a `ClothArmPose`; an analytic geometric Jacobian from the stored world axes and origins would be exact and far cheaper, but is not implemented.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 238.
