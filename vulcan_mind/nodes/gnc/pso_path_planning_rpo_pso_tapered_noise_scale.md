---
id: gnc.pso_path_planning_rpo_pso_tapered_noise_scale
label: rpo_pso_tapered_noise_scale
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_tapered_noise_scale
  lines:
  - 33
  - 33
inputs:
- id: j
  type: Int
  units: n/a
  required: true
  description: Positional argument `j`.
- id: n_waypoints
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
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
  description: Return value of `rpo_pso_tapered_noise_scale`. Returns `clamp(edge_distance
    / max((n_waypoints + 1) / 2, 1.0), 0.25, 1.0)`.
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

# rpo_pso_tapered_noise_scale

## Purpose
Reduces reseeding noise on waypoints near the path ends so culled particles keep their approach and departure geometry while their middles are shaken up.

## Design & Implementation
Returns 1.0 for a single waypoint. Otherwise it measures the waypoint's distance from the nearer end as `min(j, n+1-j)`, divides by half the span, and clamps into `[0.25, 1.0]`. A central waypoint therefore gets full noise and the outermost ones a quarter.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `j` | Int | n/a | yes | Positional argument `j`. |
| in | `n_waypoints` | Int | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_tapered_noise_scale`. Returns `clamp(edge_distance / max((n_waypoints + 1) / 2, 1.0), 0.25, 1.0)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_cull_swarm_bang|rpo_pso_cull_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:87-87`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The 0.25 floor and the linear taper are fixed, with no configuration hook; for very long polygons the central plateau at 1.0 covers most waypoints so the taper only affects the two or three at each end.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 33.
