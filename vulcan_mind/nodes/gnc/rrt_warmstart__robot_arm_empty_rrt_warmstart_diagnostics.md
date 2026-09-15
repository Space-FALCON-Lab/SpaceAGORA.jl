---
id: gnc.rrt_warmstart__robot_arm_empty_rrt_warmstart_diagnostics
label: _robot_arm_empty_rrt_warmstart_diagnostics
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_empty_rrt_warmstart_diagnostics
  lines:
  - 185
  - 185
inputs:
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
  description: Return value of `_robot_arm_empty_rrt_warmstart_diagnostics`. Returns
    `(`.
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

# _robot_arm_empty_rrt_warmstart_diagnostics

## Purpose
`_robot_arm_empty_rrt_warmstart_diagnostics` produces the diagnostics record reported when the RRT warm start is disabled or never attempted, so downstream result structures always have the same fields regardless of whether a warm start ran.

## Design & Implementation
Signature `(cfg::RobotArmHYPRConfig)`. It returns a `NamedTuple` with `enabled = cfg.rrt_warmstart_enable`, `attempted = false`, `path_found = false`, `iterations = 0`, `cost = Inf`, `raw_cost = Inf`, and `n_points = 0`. The field set matches the diagnostics constructed on the success and failure paths of `_robot_arm_rrt_connect_warmstart_path`, so `_robot_arm_rrt_warmstart_fields` can flatten either without special-casing. The function allocates nothing beyond the tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_empty_rrt_warmstart_diagnostics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:197-197`
- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:221-221`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`cost = Inf` is used to signal absence, which any consumer computing statistics over costs (means, minima) must filter out explicitly. The tuple carries no reason field, so a disabled warm start and one skipped for another reason are indistinguishable except through `enabled`. Field names are duplicated in three places in the file and must be kept in sync by hand.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 185.
