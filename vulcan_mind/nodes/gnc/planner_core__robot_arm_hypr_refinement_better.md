---
id: gnc.planner_core__robot_arm_hypr_refinement_better
label: _robot_arm_hypr_refinement_better
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: _robot_arm_hypr_refinement_better
  lines:
  - 48
  - 48
inputs:
- id: candidate
  type: Any
  units: n/a
  required: true
  description: Positional argument `candidate`.
- id: current
  type: Any
  units: n/a
  required: true
  description: Positional argument `current`.
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
  type: Bool
  units: n/a
  description: Return value of `_robot_arm_hypr_refinement_better`. Returns `true`
    or `abs_improvement > cfg.refinement_min_abs_cost_improvement &&`.
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

# _robot_arm_hypr_refinement_better

## Purpose
Lexicographic acceptance test used by `_robot_arm_hypr_post_refine_points` to decide whether a perturbed path (`candidate`) should replace the incumbent (`current`). Both are cost-component named tuples from `robot_arm_hypr_path_cost_components`. Obstacle safety dominates, then clearance while infeasible, then a thresholded improvement in total cost.

## Design & Implementation
Returns `true` immediately if `candidate.J_obs < current.J_obs`, and `false` if it is worse by more than `1e-9`. When obstacle cost is tied and the incumbent is still infeasible (`current.J_obs > 0`), a clearance gain of more than `1e-6` m is accepted. Otherwise it computes `abs_improvement = current.total - candidate.total` and `rel_improvement = abs_improvement / max(abs(current.total), 1e-12)` and requires both to exceed `cfg.refinement_min_abs_cost_improvement` and `cfg.refinement_min_rel_cost_improvement`. The function is pure and never throws.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidate` | Any | n/a | yes | Positional argument `candidate`. |
| in | `current` | Any | n/a | yes | Positional argument `current`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_robot_arm_hypr_refinement_better`. Returns `true` or `abs_improvement > cfg.refinement_min_abs_cost_improvement &&`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core__robot_arm_hypr_post_refine_points|_robot_arm_hypr_post_refine_points]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:99-99`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
With `current.total = Inf` (pruned incumbent) both improvement metrics become `Inf` or `NaN`; `Inf - Inf` yields `NaN`, so comparisons return `false` and the candidate is rejected even if better. Tolerances `1e-9` and `1e-6` are hard-coded rather than taken from `cfg`. Because a clearance gain is only rewarded while `J_obs > 0`, feasible paths cannot be refined toward larger margins.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl` line 48.
