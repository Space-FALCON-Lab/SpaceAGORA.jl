---
id: dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat
label: _rest_child_parent_quat
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _rest_child_parent_quat
  lines:
  - 160
  - 160
inputs:
- id: parent_q
  type: Any
  units: n/a
  required: true
  description: Positional argument `parent_q`.
- id: child_q
  type: Any
  units: n/a
  required: true
  description: Positional argument `child_q`.
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
  description: Return value of `_rest_child_parent_quat`. Returns `_quat_mul(_quat_conj(parent_q),
    child_q)`.
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

# _rest_child_parent_quat

## Purpose
Computes the relative orientation of a child link with respect to its parent, `conj(parent_q) * child_q`, which becomes the joint's rest quaternion.

## Design & Implementation
A one-line composition of `_quat_conj` and `_quat_mul`, both of which normalise. The result satisfies `child_q ≈ parent_q * rest`, which is exactly how `assign_coupled_cloth_robot_arm_rhs!` reconstructs `desired_child_q`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `parent_q` | Any | n/a | yes | Positional argument `parent_q`. |
| in | `child_q` | Any | n/a | yes | Positional argument `child_q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rest_child_parent_quat`. Returns `_quat_mul(_quat_conj(parent_q), child_q)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:301-301`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:172-172`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- `callees` → [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- `callees` → [[vehicle.robotics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
<!-- vulcan:connections:end -->

## Limitations
Four normalisations per call (two inside `_quat_mul`, one in `_quat_conj`, one on output). Degenerate inputs silently become identity. No sign canonicalisation is applied, so `rest` and `-rest` may both appear across time samples.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 160.
