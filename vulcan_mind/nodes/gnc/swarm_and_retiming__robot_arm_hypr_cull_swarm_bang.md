---
id: gnc.swarm_and_retiming__robot_arm_hypr_cull_swarm_bang
label: _robot_arm_hypr_cull_swarm!
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_cull_swarm!
  lines:
  - 21
  - 21
inputs:
- id: positions
  type: Any
  units: n/a
  required: true
  description: Positional argument `positions`.
- id: velocities
  type: Any
  units: n/a
  required: true
  description: Positional argument `velocities`.
- id: pbest
  type: Any
  units: n/a
  required: true
  description: Positional argument `pbest`.
- id: pbest_cost
  type: Any
  units: n/a
  required: true
  description: Positional argument `pbest_cost`.
- id: lo_rep
  type: Any
  units: n/a
  required: true
  description: Positional argument `lo_rep`.
- id: hi_rep
  type: Any
  units: n/a
  required: true
  description: Positional argument `hi_rep`.
- id: gbest
  type: Any
  units: n/a
  required: true
  description: Positional argument `gbest`.
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
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
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
  type: Int
  units: n/a
  description: Return value of `_robot_arm_hypr_cull_swarm!`; mutates `positions`
    in place. Returns `0` or `n_replace`.
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

# _robot_arm_hypr_cull_swarm!

## Purpose
Replaces the worst-performing fraction of HYPR particles with jittered copies of the global best, restarting exploration near the incumbent solution once culling is enabled.

## Design & Implementation
Returns 0 without mutation when `cfg.cull_enable` is false, `iter < cfg.cull_start_iter`, `cfg.cull_fraction_max <= 0`, or `minimum(pbest_cost)` is non-finite. Otherwise `n_replace = min(n_particles - 1, floor(Int, cull_fraction_max * n_particles))`; the `n_replace` particles with the highest `pbest_cost` are chosen by `sortperm(...; rev=true)`. For each dimension `d`, the new position is `clamp(gbest[d] + cull_noise_scale * (hi_rep[d] - lo_rep[d]) * randn(rng), lo_rep[d], hi_rep[d])`; velocity is zeroed, `pbest` is set to the new position, and `pbest_cost` is set to `Inf`. Mutates `positions`, `velocities`, `pbest`, and `pbest_cost` in place and returns `n_replace`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `positions` | Any | n/a | yes | Positional argument `positions`. |
| in | `velocities` | Any | n/a | yes | Positional argument `velocities`. |
| in | `pbest` | Any | n/a | yes | Positional argument `pbest`. |
| in | `pbest_cost` | Any | n/a | yes | Positional argument `pbest_cost`. |
| in | `lo_rep` | Any | n/a | yes | Positional argument `lo_rep`. |
| in | `hi_rep` | Any | n/a | yes | Positional argument `hi_rep`. |
| in | `gbest` | Any | n/a | yes | Positional argument `gbest`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `iter` | Int | n/a | yes | Positional argument `iter`. |
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_robot_arm_hypr_cull_swarm!`; mutates `positions` in place. Returns `0` or `n_replace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:306-306`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:306-306`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Always keeps at least one particle uncullled but never protects the global-best particle explicitly, relying on it having the lowest cost. `sortperm` allocates every call. Setting `pbest_cost` to `Inf` means culled particles cannot contribute to early-stopping statistics until re-evaluated. Bounds arrays `lo_rep`/`hi_rep` are assumed to have the same length as the position dimension.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 21.
