---
id: gnc.config__validate_robot_arm_hypr_config
label: _validate_robot_arm_hypr_config
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/config.jl
  symbol: _validate_robot_arm_hypr_config
  lines:
  - 93
  - 93
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
  description: Return value of `_validate_robot_arm_hypr_config`. Returns `cfg`.
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

# _validate_robot_arm_hypr_config

## Purpose
`_validate_robot_arm_hypr_config(cfg::RobotArmHYPRConfig)` is the single fail-fast gate applied to a robot-arm HYPR configuration before any planning work begins. It converts silently degenerate settings — zero particles, a one-sample path, an unknown curve type, a shrink factor outside the open unit interval — into an immediate, named `ArgumentError` instead of a downstream division by zero or an empty swarm.

## Design & Implementation
The body is a flat sequence of short-circuit assertions of the form `predicate || throw(ArgumentError("<field> must be ..."))`, grouped by subsystem: swarm sizing (`n_waypoints >= 0`, `n_particles > 0`, `n_iters > 0`, `n_samples >= 2`), curve and clearance (`curve_type in (:bezier, :polyline)`, `safe_distance_m >= 0`), the RRT warm-start block (iteration count, positive `rrt_warmstart_step_size_rad`, `goal_sample_rate` inside `[0,1]`, positive collision sampling step, positive connect step cap, non-negative shortcut iterations and runtime limit), the retiming block (positive joint velocity and acceleration caps, non-negative reaction time and gain, `retime_min_scale >= 1.0`, `retime_max_scale >= retime_min_scale`, positive base force and torque caps, positive wrench model gain, `retime_base_wrench_margin >= 1.0`), the refinement block (non-negative rounds, positive `refinement_step_fraction`, `0 < refinement_shrink < 1`, non-negative absolute and relative improvement thresholds) and the cloth model (positive `retime_cloth_dt_s`, non-negative stiffness and damping coefficients). On success it returns `cfg` unchanged, so it can be used inline where the configuration is consumed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_validate_robot_arm_hypr_config`. Returns `cfg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:184-184`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation is per-field and local; no cross-field consistency is checked beyond the two scale bounds, so an `early_stopping_min_iters` larger than `n_iters` or a `cull_start_iter` beyond `n_iters` passes silently and merely disables the feature. Several documented fields are never examined at all — the PSO weights `w_len`, `w_smooth`, `w_obs`, `w_inertia`, `c1`, `c2`, the schedule block, and every `cull_*` and `early_stopping_*` threshold — so negative weights or an inverted schedule reach the optimiser unchecked. `NaN` defeats every comparison used here and therefore passes every assertion.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/config.jl` line 93.
