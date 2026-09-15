---
id: vehicle.robotics_cloth_fk_state
label: cloth_fk_state
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: cloth_fk_state
  lines:
  - 212
  - 212
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
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
- id: dq
  type: Any
  units: n/a
  required: false
  description: Keyword argument `dq` (default `zeros(length(model.joints))`).
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
  type: ClothArmState
  units: n/a
  description: Return value of `cloth_fk_state`.
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

# cloth_fk_state

## Purpose
Runs forward kinematics and packages the result with joint rates into a `ClothArmState` for the dynamics layer.

## Design & Implementation
Validates both `q` and `dq` — the latter defaulting to zeros — through `_validate_joint_vector`, calls `cloth_fk`, and constructs the state with link linear and angular velocity vectors zero-filled to the link count. The zeros are placeholders the multibody code overwrites.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `dq` | Any | n/a | no | Keyword argument `dq` (default `zeros(length(model.joints))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmState | n/a | — | Return value of `cloth_fk_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_reference_state|cloth_reference_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:31-31`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/vehicle/robotics/robotics.jl:220-220`
- `callees` → [[vehicle.robotics__validate_joint_vector|_validate_joint_vector]] · `callers` · call · `src/vehicle/robotics/robotics.jl:218-218`
- `callees` → [[vehicle.robotics_clotharmstate|ClothArmState]] · `callers` · call · `src/vehicle/robotics/robotics.jl:221-221`
<!-- vulcan:connections:end -->

## Limitations
The returned link velocities are always zero regardless of `dq`; the function does not compute the velocity kinematics, so a consumer expecting them from this call alone gets a resting arm.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 212.
