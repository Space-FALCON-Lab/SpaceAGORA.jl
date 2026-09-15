---
id: gnc.swarm_and_retiming__robot_arm_hypr_iteration_weights
label: _robot_arm_hypr_iteration_weights
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_iteration_weights
  lines:
  - 2
  - 2
inputs:
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: iter
  type: Int
  units: n/a
  required: true
  description: Positional argument `iter`.
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
  description: Return value of `_robot_arm_hypr_iteration_weights`. Returns `hypr_iteration_weights(`.
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

# _robot_arm_hypr_iteration_weights

## Purpose
Provides the inertia and cognitive/social coefficients for a given HYPR PSO iteration, honouring the optional annealing schedule in the robot-arm configuration.

## Design & Implementation
Thin adapter that unpacks thirteen fields from `cfg::RobotArmHYPRConfig` and forwards them positionally to the shared `hypr_iteration_weights`: `schedule_enable`, `n_iters`, `iter`, `w_inertia`, `c1`, `c2`, `schedule_transition_fraction`, `schedule_w_min`, `schedule_w_end_fraction`, `schedule_c1_end_fraction`, `schedule_c2_end_fraction`, `schedule_c_min`, `schedule_c_max`. The returned tuple of `(w, c1, c2)` is consumed once per iteration by the swarm update loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `iter` | Int | n/a | yes | Positional argument `iter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_iteration_weights`. Returns `hypr_iteration_weights(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:319-319`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:319-319`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_iteration_weights|hypr_iteration_weights]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
Argument order is positional and must stay in sync with `hypr_iteration_weights`; a reordering in the shared function would silently produce wrong weights. No validation of `iter` against `cfg.n_iters` occurs here.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 2.
