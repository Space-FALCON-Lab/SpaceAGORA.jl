---
id: dynamics.cloth_multibody_compliant_state_parts
label: compliant_state_parts
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: compliant_state_parts
  lines:
  - 458
  - 458
inputs:
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: Return value of `compliant_state_parts`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# compliant_state_parts

## Purpose
Unpacks one body's block of the flat state into a named tuple of position, unit quaternion, velocity and angular rate.

## Design & Implementation
Computes the block offset `13(i-1)` and slices four static vectors, renormalising the quaternion on the way out so consumers always see a unit orientation even mid-integration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compliant_state_parts`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__parent_kinematics|_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:496-496`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_end_effector|cloth_robot_arm_end_effector]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:430-430`
- [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:615-615`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:462-462`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:462-462`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:462-462`
<!-- vulcan:connections:end -->

## Limitations
Slicing an `AbstractVector` with ranges allocates before the `SVector` conversion in some Julia versions; called several times per body per derivative.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 458.
