---
id: gnc.rrt_warmstart__robot_arm_rrt_warmstart_fields
label: _robot_arm_rrt_warmstart_fields
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_warmstart_fields
  lines:
  - 198
  - 198
inputs:
- id: warmstart
  type: Any
  units: n/a
  required: true
  description: Positional argument `warmstart`.
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
  description: Return value of `_robot_arm_rrt_warmstart_fields`. Returns `(`.
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

# _robot_arm_rrt_warmstart_fields

## Purpose
`_robot_arm_rrt_warmstart_fields` renames the warm-start diagnostics into `rrt_warmstart_`-prefixed keys and attaches the original tuple, producing the flat set of fields that is merged into the `RobotArmHYPRResult` diagnostics so that warm-start metrics sit alongside optimiser metrics without name clashes.

## Design & Implementation
Signature `(warmstart)`, accepting any object with the seven diagnostic properties. It returns a `NamedTuple` mapping `rrt_warmstart_enabled`, `rrt_warmstart_attempted`, `rrt_warmstart_path_found`, `rrt_warmstart_iterations`, `rrt_warmstart_cost`, `rrt_warmstart_raw_cost` and `rrt_warmstart_n_points` to the corresponding `warmstart.*` values, plus `rrt_warmstart = warmstart` holding the unflattened record. Pure and allocation-light; used once per plan.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `warmstart` | Any | n/a | yes | Positional argument `warmstart`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_warmstart_fields`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:351-351`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:217-217`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The input is untyped, so a tuple missing one of the seven properties fails with an `ErrorException` on field access at call time rather than a clear validation error. Keeping the nested `rrt_warmstart` alongside the flattened keys duplicates data in the result. Any new diagnostic added to the warm-start tuple must also be added here or it is silently dropped from the flat view.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 198.
