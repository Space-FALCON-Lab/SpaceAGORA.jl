---
id: vehicle.robotics_cloth_end_effector_pose
label: cloth_end_effector_pose
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: cloth_end_effector_pose
  lines:
  - 231
  - 231
inputs:
- id: state
  type: ClothArmState
  units: n/a
  required: true
  description: Positional argument `state`.
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
  description: Return value of `cloth_end_effector_pose`.
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

# cloth_end_effector_pose

## Purpose
Convenience accessor returning the end-effector position and quaternion as a pair from either a pose or a state.

## Design & Implementation
Two one-line methods: the `ClothArmState` method reaches through `state.pose`, the `ClothArmPose` method reads the fields directly. Both return a tuple of the static position and the static scalar-last quaternion, so a caller can destructure without knowing which type it was handed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | ClothArmState | n/a | yes | Positional argument `state`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_end_effector_pose`. |
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
Only the last link's tip is exposed; a tool frame offset from the tip must be applied by the caller using the returned quaternion.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 231.
