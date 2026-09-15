---
id: gnc.trajectory_optimizers_rpo_stomp_waypoint_state_cost
label: rpo_stomp_waypoint_state_cost
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_stomp_waypoint_state_cost
  lines:
  - 336
  - 336
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: internal_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `internal_idx`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: safe_distance_m
  type: Any
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
- id: w_obs
  type: Any
  units: n/a
  required: true
  description: Keyword argument `w_obs`.
- id: w_len
  type: Any
  units: n/a
  required: true
  description: Keyword argument `w_len`.
- id: w_smooth
  type: Any
  units: n/a
  required: true
  description: Keyword argument `w_smooth`.
- id: obstacle_margin_m
  type: Any
  units: n/a
  required: true
  description: Keyword argument `obstacle_margin_m`.
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
  description: Return value of `rpo_stomp_waypoint_state_cost`. Returns `w_obs * obs
    + w_len * local_len + w_smooth * dot(d2, d2)`.
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

# rpo_stomp_waypoint_state_cost

## Purpose
Computes the local per-waypoint cost STOMP uses to weight rollouts independently at each internal waypoint: obstacle potential at the point, the length of its two adjacent segments, and its squared second difference.

## Theory & Math
$$S_i = w_{obs}\,c(d_i) + w_{len}\left(\|p_i - p_{i-1}\| + \|p_{i+1} - p_i\|\right) + w_{smooth}\,\|p_{i-1} - 2p_i + p_{i+1}\|^2$$ where $c(\cdot)$ is the CHOMP potential and $p_i$ is the waypoint at full-trajectory index $i$.

## Design & Implementation
Signature `rpo_stomp_waypoint_state_cost(points, internal_idx, geometry; safe_distance_m, w_obs, w_len, w_smooth, obstacle_margin_m)`. Maps `point_idx = internal_idx + 1` into the full trajectory (skipping the start column), queries `rpo_clearance_distance_to_station`, evaluates `rpo_chomp_obstacle_potential`, sums `norm(p - prev) + norm(next - p)`, and forms `d2 = prev - 2p + next`. Returns `w_obs * obs + w_len * local_len + w_smooth * dot(d2, d2)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `internal_idx` | Integer | n/a | yes | Positional argument `internal_idx`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Any | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `w_obs` | Any | n/a | yes | Keyword argument `w_obs`. |
| in | `w_len` | Any | n/a | yes | Keyword argument `w_len`. |
| in | `w_smooth` | Any | n/a | yes | Keyword argument `w_smooth`. |
| in | `obstacle_margin_m` | Any | n/a | yes | Keyword argument `obstacle_margin_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_stomp_waypoint_state_cost`. Returns `w_obs * obs + w_len * local_len + w_smooth * dot(d2, d2)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:425-425`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:348-348`
- `callees` → [[gnc.trajectory_optimizers_rpo_chomp_obstacle_potential|rpo_chomp_obstacle_potential]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:349-349`
<!-- vulcan:connections:end -->

## Limitations
Unlike the global objective, the local terms are un-normalised (raw metres), so `w_len` and `w_obs` act on different scales here than in `rpo_trajectory_soft_objective`. Only the waypoint itself is checked for clearance, not the segments between waypoints, so a rollout can thread through the station between samples and still score low. No bounds check on `internal_idx`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 336.
