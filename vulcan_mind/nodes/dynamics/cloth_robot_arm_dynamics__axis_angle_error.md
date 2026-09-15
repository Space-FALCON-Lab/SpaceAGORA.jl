---
id: dynamics.cloth_robot_arm_dynamics__axis_angle_error
label: _axis_angle_error
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _axis_angle_error
  lines:
  - 135
  - 135
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
Converts an orientation-error quaternion into a 3-vector rotation error (axis times angle, radians) that the rotational spring and PD actuator can act on.

## Theory & Math
For unit $q_{err} = (\mathbf{v}, w)$ with $w \ge 0$: $\theta = 2\,\operatorname{atan2}(\lVert\mathbf{v}\rVert, w)$ and $\boldsymbol\phi = \theta\,\mathbf{v}/\lVert\mathbf{v}\rVert$, reducing to $\boldsymbol\phi \approx 2\mathbf{v}$ as $\theta \to 0$.

## Design & Implementation
Normalises with `_unit_quat`, flips sign when the scalar part `q[4] < 0` so the shortest rotation is chosen, extracts the vector part `v`, and computes `nv = norm(v)`. When `nv <= 1e-12` it returns the small-angle approximation `2v`; otherwise returns `(θ / nv) * v` with `θ = 2 * atan(nv, q[4])`.

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

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:136-136`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:136-136`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:136-136`
<!-- vulcan:connections:end -->

## Limitations
The output is in whatever frame `q_err` was composed in (world frame in the RHS), and callers must rotate consistently. Errors of exactly 180 degrees have `w = 0` and an arbitrary sign choice. The `1e-12` switch threshold is hard-coded.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 135.
