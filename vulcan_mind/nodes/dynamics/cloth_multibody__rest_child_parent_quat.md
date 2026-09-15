---
id: dynamics.cloth_multibody__rest_child_parent_quat
label: _rest_child_parent_quat
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _rest_child_parent_quat
  lines:
  - 276
  - 276
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
Computes the relative orientation of a child body with respect to its parent, used as the unloaded orientation of a joint's rotational spring when none is given explicitly.

## Design & Implementation
Returns `conj(parent_q) ⊗ child_q` through the normalising multiply. Called by `build_compliant_topology` from the nodes' initial quaternions, so a topology built in its rest configuration has zero rotational preload.

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
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- `callees` → [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- `callees` → [[vehicle.robotics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
<!-- vulcan:connections:end -->

## Limitations
If the initial configuration is not the intended rest configuration the springs are preloaded from the first step, and the only remedy is to pass explicit rest quaternions on the edges.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 276.
