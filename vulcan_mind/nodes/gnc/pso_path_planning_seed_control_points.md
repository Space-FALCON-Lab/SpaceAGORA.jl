---
id: gnc.pso_path_planning_seed_control_points
label: seed_control_points
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: seed_control_points
  lines:
  - 290
  - 290
inputs:
- id: seed_points
  type: Any
  units: n/a
  required: true
  description: Positional argument `seed_points`.
- id: n_waypoints
  type: Any
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
- id: local_cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `local_cfg`.
- id: seed_curve_type
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `seed_curve_type` (default `local_cfg.curve_type`).
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
  description: Return value of `seed_control_points`. Returns `points` or `seeded`.
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

# seed_control_points

## Purpose
Produces the initial control polygon the first particle is placed on, from a warm-start path if one exists or the straight start-goal chord otherwise.

## Design & Implementation
A closure. With no seed it distributes `n_waypoints` interior points evenly along the chord between the captured `start` and `goal`. With a seed and a Bezier curve type it samples the seed at `sample_ds_m` — as a Bezier or polyline depending on `seed_curve_type` — and fits a fixed-endpoint Bezier with `n_waypoints + 2` controls; for polyline curve types it resamples the seed to that many points. The first and last columns are overwritten with the exact start and goal.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `seed_points` | Any | n/a | yes | Positional argument `seed_points`. |
| in | `n_waypoints` | Any | n/a | yes | Positional argument `n_waypoints`. |
| in | `local_cfg` | RPOPSOConfig | n/a | yes | Positional argument `local_cfg`. |
| in | `seed_curve_type` | Symbol | n/a | no | Keyword argument `seed_curve_type` (default `local_cfg.curve_type`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `seed_control_points`. Returns `points` or `seeded`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:342-342`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:290-290`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_resample_polyline_points|rpo_resample_polyline_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:307-307`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:303-303`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:304-304`
- `callees` → [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:305-305`
<!-- vulcan:connections:end -->

## Limitations
Fitting a low-degree Bezier to a jagged RRT polyline smooths away detours the RRT took for clearance, so the seed can be infeasible even though the RRT path was not.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 290.
