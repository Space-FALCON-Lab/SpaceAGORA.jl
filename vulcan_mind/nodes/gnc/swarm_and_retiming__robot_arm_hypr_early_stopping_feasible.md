---
id: gnc.swarm_and_retiming__robot_arm_hypr_early_stopping_feasible
label: _robot_arm_hypr_early_stopping_feasible
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_early_stopping_feasible
  lines:
  - 64
  - 64
inputs:
- id: components
  type: Any
  units: n/a
  required: true
  description: Positional argument `components`.
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_robot_arm_hypr_early_stopping_feasible`. Returns
    `getproperty(components, :J_obs) <= 1.0e-9`.
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

# _robot_arm_hypr_early_stopping_feasible

## Purpose
Decides whether the current best HYPR solution is feasible enough for early stopping to be allowed, based on the obstacle-penalty component of the cost.

## Design & Implementation
If `cfg.early_stopping_require_feasible` is false the function returns `true` unconditionally. Otherwise it reads `getproperty(components, :J_obs)` from the cost-components record and returns `J_obs <= 1.0e-9`, treating any obstacle penalty above that tolerance as infeasible.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `components` | Any | n/a | yes | Positional argument `components`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_early_stopping_feasible`. Returns `getproperty(components, :J_obs) <= 1.0e-9`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:291-291`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:291-291`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the obstacle term is inspected; joint-limit or wrench violations encoded in other components do not block early stopping. The tolerance `1.0e-9` is hard-coded and absolute, so cost scaling changes the meaning of feasibility. `components` must expose a `J_obs` property or `getproperty` throws.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 64.
