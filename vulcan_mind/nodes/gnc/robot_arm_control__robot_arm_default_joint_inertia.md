---
id: gnc.robot_arm_control__robot_arm_default_joint_inertia
label: _robot_arm_default_joint_inertia
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: _robot_arm_default_joint_inertia
  lines:
  - 37
  - 37
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
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
  description: Return value of `_robot_arm_default_joint_inertia`. Returns `inertia`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _robot_arm_default_joint_inertia

## Purpose
Provides a rough per-joint effective inertia for the MPC prediction model when the user does not supply `joint_inertia_kg_m2`, by treating every link distal to a joint as a point mass at its centre-of-mass distance from that joint along the stretched-out chain.

## Theory & Math
For joint $i$ with distal links $j = i,\dots,n$: $J_i = \max\!\left(\sum_{j=i}^{n} m_j \max(d_j, r_j)^2,\ 10^{-6}\right)$, where $m_j$ is `mass_kg`, $r_j$ is `radius_m`, and $d_j = \sum_{k=i}^{j-1}\|\mathbf{v}_k\| + \|\mathbf{c}_j\|$ is the centre-of-mass distance along the chain built from the link vectors $\mathbf{v}_k$ (`vector_parent`) and the COM offset $\mathbf{c}_j$ (`com_offset_parent`).

## Design & Implementation
Signature `_robot_arm_default_joint_inertia(plan::RobotArmPlan)`. With `links = plan.model.links` and `n = length(links)`, the output starts as `fill(1.0e-4, n)`. For joint `i` it walks `j = i:n`, accumulating `distance_to_joint` by `norm(link.vector_parent)` and adding `link.mass_kg * max(com_distance, link.radius_m)^2` where `com_distance = distance_to_joint + norm(link.com_offset_parent)`. The floor `max(total, 1.0e-6)` keeps inertias strictly positive for the `Diagonal(1 ./ inertia)` inversion in `init_robot_arm_joint_mpc`. Returns a `Vector{Float64}` in kg m^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_default_joint_inertia`. Returns `inertia`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:72-72`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The estimate assumes a fully extended, collinear arm and point-mass links, ignoring joint axes, link orientation and the link's own rotational inertia, so it can be off by a large factor for folded configurations; the MPC then tracks with a mismatched plant model. The magnitudes `1.0e-4` and `1.0e-6` are arbitrary floors with no physical basis. The nested loop is O(n^2), negligible for small arms but recomputed on every `init_robot_arm_joint_mpc` call.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 37.
