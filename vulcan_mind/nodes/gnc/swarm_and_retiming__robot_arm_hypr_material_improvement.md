---
id: gnc.swarm_and_retiming__robot_arm_hypr_material_improvement
label: _robot_arm_hypr_material_improvement
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_material_improvement
  lines:
  - 54
  - 54
inputs:
- id: new_cost
  type: Real
  units: n/a
  required: true
  description: Positional argument `new_cost`.
- id: reference_cost
  type: Real
  units: n/a
  required: true
  description: Positional argument `reference_cost`.
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
  description: Return value of `_robot_arm_hypr_material_improvement`. Returns `hypr_material_improvement(`.
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

# _robot_arm_hypr_material_improvement

## Purpose
Tests whether a new HYPR cost value is a large enough improvement over a reference to reset the early-stopping stall counter.

## Design & Implementation
Adapter around the shared `hypr_material_improvement(new_cost, reference_cost, min_abs, min_rel)` supplying `cfg.early_stopping_min_abs_improvement` and `cfg.early_stopping_min_rel_improvement` from `cfg::RobotArmHYPRConfig`. Both inputs are `Real`; the return is the boolean produced by the shared routine, which requires the improvement to clear both the absolute and relative thresholds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `new_cost` | Real | n/a | yes | Positional argument `new_cost`. |
| in | `reference_cost` | Real | n/a | yes | Positional argument `reference_cost`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_material_improvement`. Returns `hypr_material_improvement(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:292-292`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:292-292`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_material_improvement|hypr_material_improvement]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
The semantics (whether both thresholds must be met or either suffices) live entirely in `hypr_material_improvement`; nothing here documents or enforces them. Non-finite costs are passed straight through.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 54.
