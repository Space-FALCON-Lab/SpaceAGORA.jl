---
id: gnc.pso_path_planning_reset_swarm_bang
label: reset_swarm!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: reset_swarm!
  lines:
  - 315
  - 315
inputs:
- id: new_n_waypoints
  type: Any
  units: n/a
  required: true
  description: Positional argument `new_n_waypoints`.
- id: new_search_margin
  type: Any
  units: n/a
  required: true
  description: Positional argument `new_search_margin`.
- id: seed_points
  type: Any
  units: n/a
  required: false
  description: Keyword argument `seed_points` (default `nothing`).
- id: seed_curve_type
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `seed_curve_type` (default `cfg.curve_type`).
- id: use_warmstart_bounds
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `use_warmstart_bounds` (default `false`).
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
  type: Nothing
  units: n/a
  description: Return value of `reset_swarm!`; mutates `new_n_waypoints` in place.
    Returns `nothing`.
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

# reset_swarm!

## Purpose
Reinitialises every particle array for a new waypoint count and search margin, placing particle one on the seed and the rest in a Gaussian cloud around it.

## Design & Implementation
A closure that updates the captured `current_n_waypoints` and `current_search_margin`, computes bounds either from the warm-start path or from start and goal, and reallocates `positions`, `velocities`, `pbest`, `pbest_cost`, `pbest_obs`, `stagnation_count` and `gbest`. It flattens the seeded polygon's interior points into `base`, then sets particle one to `base` clamped, and every other particle to `base` plus `spread_scale` times the per-dimension span times a standard normal, clamped. Velocities start at a tenth of a uniform random fraction of the span.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `new_n_waypoints` | Any | n/a | yes | Positional argument `new_n_waypoints`. |
| in | `new_search_margin` | Any | n/a | yes | Positional argument `new_search_margin`. |
| in | `seed_points` | Any | n/a | no | Keyword argument `seed_points` (default `nothing`). |
| in | `seed_curve_type` | Symbol | n/a | no | Keyword argument `seed_curve_type` (default `cfg.curve_type`). |
| in | `use_warmstart_bounds` | Bool | n/a | no | Keyword argument `use_warmstart_bounds` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `reset_swarm!`; mutates `new_n_waypoints` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:498-498`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:315-315`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:323-323`
- `callees` → [[gnc.pso_helpers_rpo_pso_warmstart_bounds|rpo_pso_warmstart_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:326-326`
- `callees` → [[gnc.pso_path_planning_cfg_for_current|cfg_for_current]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:324-324`
- `callees` → [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:342-342`
- `callees` → [[gncy.pso_helpers_rpo_pso_bounds|rpo_pso_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:327-327`
<!-- vulcan:connections:end -->

## Limitations
It rebinds captured variables rather than mutating arrays in place, which works because every other closure reads the same captured bindings, but it means holding a reference to the old `positions` array after a reset observes stale data.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 315.
