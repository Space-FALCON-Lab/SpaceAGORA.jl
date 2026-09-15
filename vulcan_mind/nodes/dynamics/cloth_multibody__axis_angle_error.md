---
id: dynamics.cloth_multibody__axis_angle_error
label: _axis_angle_error
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _axis_angle_error
  lines:
  - 165
  - 165
inputs:
- id: q_err_raw
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_err_raw`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_axis_angle_error`.
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

# _axis_angle_error

## Purpose
Converts an orientation error quaternion into the rotation vector the rotational spring and PD actuator act on.

## Theory & Math
$$
\theta = 2\,\operatorname{atan2}(\|\vec{v}\|, w),\qquad \vec{\phi} = \frac{\theta}{\|\vec{v}\|}\,\vec{v}
$$

## Design & Implementation
Normalises the error, flips its sign if the scalar part is negative so the shortest rotation is chosen, extracts the vector part `v`, and returns `2v` for tiny norms or `(θ / |v|) v` with `θ = 2 atan(|v|, w)` otherwise. The `atan` form is accurate across the whole range where `acos` would lose precision near zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_err_raw` | Any | n/a | yes | Positional argument `q_err_raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_axis_angle_error`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:549-549`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:394-394`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:166-166`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:166-166`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
At exactly 180 degrees the axis is ambiguous and the sign flip picks one arbitrarily, so a joint starting fully inverted can spring either way.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 165.
