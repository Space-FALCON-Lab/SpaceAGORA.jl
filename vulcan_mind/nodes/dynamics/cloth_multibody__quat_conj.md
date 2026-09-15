---
id: dynamics.cloth_multibody__quat_conj
label: _quat_conj
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _quat_conj
  lines:
  - 122
  - 122
inputs:
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  type: SVector{4,
  units: n/a
  description: Return value of `_quat_conj`.
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

# _quat_conj

## Purpose
Returns the conjugate, hence inverse, of a unit quaternion.

## Design & Implementation
Normalises the input then negates the vector part, keeping the scalar. `@inline`. Used to form relative orientations between parent and child.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `_quat_conj`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__rest_child_parent_quat|_rest_child_parent_quat]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:548-548`
- [[dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat|_rest_child_parent_quat]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:393-393`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:123-123`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:123-123`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:123-123`
<!-- vulcan:connections:end -->

## Limitations
Normalising first means the result is always the inverse of the normalised input, which differs from the true inverse of a non-unit quaternion.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 122.
