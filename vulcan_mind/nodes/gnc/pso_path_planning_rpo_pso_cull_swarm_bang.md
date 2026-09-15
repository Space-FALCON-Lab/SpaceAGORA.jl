---
id: gnc.pso_path_planning_rpo_pso_cull_swarm_bang
label: rpo_pso_cull_swarm!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_cull_swarm!
  lines:
  - 40
  - 40
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
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
- id: n_waypoints
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
- id: cfg
  type: RPOPSOConfig
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
  description: Return value of `rpo_pso_cull_swarm!`; mutates `positions` in place.
    Returns `0` or `n_replace`.
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

# rpo_pso_cull_swarm!

## Purpose
Replaces the worst-scoring fraction of the swarm with fresh particles biased toward a straightened version of the global best, restoring diversity late in a run.

## Design & Implementation
Returns zero unless culling is enabled, the iteration has reached `cull_start_iter`, the fraction is positive and a finite best cost exists. It selects the `n_replace` worst particles by `pbest_cost`, capped at `n_particles - 1` so at least one elite survives. For non-Bezier or waypoint-free problems it jitters each around `gbest` uniformly within `cull_noise_scale` of the box span and zeroes velocity. For Bezier paths it converts `gbest` to a path and, per waypoint, forms a target as the current point plus a 0.6-weighted pull toward the start-goal chord and a 0.4-weighted pull toward the local neighbour chord, adds tapered Gaussian noise, clamps into bounds, and sets velocity to `cull_arc_velocity_scale` times the displacement. Every replaced particle's `pbest_cost` is reset to `Inf`.

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
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `n_waypoints` | Int | n/a | yes | Positional argument `n_waypoints`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `iter` | Int | n/a | yes | Positional argument `iter`. |
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `rpo_pso_cull_swarm!`; mutates `positions` in place. Returns `0` or `n_replace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:546-546`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:546-546`

**Downstream**

- `callees` → [[gnc.pso_helpers_rpo_position_to_path|rpo_position_to_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:77-77`
- `callees` → [[gnc.pso_path_planning_rpo_pso_project_to_segment|rpo_pso_project_to_segment]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:84-84`
- `callees` → [[gnc.pso_path_planning_rpo_pso_tapered_noise_scale|rpo_pso_tapered_noise_scale]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:87-87`
<!-- vulcan:connections:end -->

## Limitations
Mutates `positions`, `velocities`, `pbest` and `pbest_cost` in place. The 0.6 and 0.4 blend weights are literals; and because `pbest_cost` is set to `Inf`, a replaced particle's previous best is forgotten even if it was better than its new location.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 40.
