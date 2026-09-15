---
id: gnc.swarm_and_retiming__robot_arm_hypr_cloth_state_for_reaction
label: _robot_arm_hypr_cloth_state_for_reaction
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_cloth_state_for_reaction
  lines:
  - 215
  - 215
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: parts_fn
  type: Any
  units: n/a
  required: true
  description: Positional argument `parts_fn`.
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
  description: Return value of `_robot_arm_hypr_cloth_state_for_reaction`. Returns
    `state`.
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

# _robot_arm_hypr_cloth_state_for_reaction

## Purpose
Packs a cloth-multibody state vector into the NamedTuple layout that `assign_coupled_cloth_robot_arm_rhs!` expects, with a placeholder spacecraft base so only the arm links carry real data.

## Design & Implementation
Takes `plan::RobotArmPlan`, the raw simulation state `x`, and `parts_fn` (bound to `ClothMultibody.compliant_state_parts`). With `n = length(plan.model.links)` it allocates `pos`, `vel`, `ω` as `zeros(3)`, `q = [0,0,0,1]` (identity quaternion, scalar-last), `mass = 1.0`, `heat_loads = zeros(1)`, and the per-link matrices `arm_r::3xn`, `arm_q::4xn`, `arm_v::3xn`, `arm_ω::3xn`. For each link `i` it calls `parts_fn(x, i)` and copies `part.r`, `part.q`, `part.v`, `part.ω` into column `i` under `@inbounds`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `parts_fn` | Any | n/a | yes | Positional argument `parts_fn`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_cloth_state_for_reaction`. Returns `state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:284-284`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The base spacecraft state is fictitious (rest at origin, unit mass), so the returned tuple is only valid for extracting arm-reaction wrenches, not for propagating the coupled system. Allocation of ten arrays per call makes it unsuitable for hot loops. `parts_fn` must return fields named exactly `r`, `q`, `v`, `ω`; no validation is performed.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 215.
