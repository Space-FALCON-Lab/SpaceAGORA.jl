---
id: vehicle.robotics_clotharmbasepose
label: ClothArmBasePose
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmBasePose
  lines:
  - 14
  - 14
inputs:
- id: position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `position`.
- id: quaternion
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `quaternion`.
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
  type: ClothArmBasePose
  units: n/a
  description: Constructed `ClothArmBasePose`.
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

# ClothArmBasePose

## Purpose
Where the arm is mounted: the world-frame position and orientation of the chain's root.

## Design & Implementation
An immutable struct of a static position and a static scalar-last quaternion. Two outer constructors accept a position as any real vector or 3-tuple, defaulting the quaternion to identity, or a position and quaternion pair, with the quaternion passed through `_unit_quat` so a slightly non-unit input is normalised on the way in.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `position` | SVector{3, Float64} | n/a | yes | Field `position`. |
| in | `quaternion` | SVector{4, Float64} | n/a | yes | Field `quaternion`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmBasePose | n/a | — | Constructed `ClothArmBasePose`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:209-209`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:26-26`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:26-26`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:26-26`
<!-- vulcan:connections:end -->

## Limitations
A zero or non-finite quaternion is silently replaced by identity rather than rejected, so a corrupted orientation input produces an upright arm with no diagnostic.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 14.
