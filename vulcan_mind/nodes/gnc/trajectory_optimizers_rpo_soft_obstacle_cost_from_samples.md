---
id: gnc.trajectory_optimizers_rpo_soft_obstacle_cost_from_samples
label: rpo_soft_obstacle_cost_from_samples
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_soft_obstacle_cost_from_samples
  lines:
  - 169
  - 169
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `rpo_soft_obstacle_cost_from_samples`. Returns `acc
    / n_samples`.
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

# rpo_soft_obstacle_cost_from_samples

## Purpose
Averages the CHOMP soft obstacle potential over all densely sampled points of a candidate path, giving the `J_obs_soft` term used by the trajectory optimizer objective in place of HYPR's hard collision count.

## Design & Implementation
Signature `rpo_soft_obstacle_cost_from_samples(samples, geometry; safe_distance_m, obstacle_margin_m)`. Returns `0.0` for zero columns. For each column `i` of `samples` it calls `rpo_clearance_distance_to_station(samples[:, i], geometry)` and accumulates `rpo_chomp_obstacle_potential(clearance, safe_distance_m, obstacle_margin_m)`, finally dividing by `n_samples` so the result is a mean potential in metres rather than a sum that grows with sampling density.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Any | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `obstacle_margin_m` | Any | n/a | yes | Keyword argument `obstacle_margin_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_soft_obstacle_cost_from_samples`. Returns `acc / n_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:194-194`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:175-175`
- `callees` → [[gnc.trajectory_optimizers_rpo_chomp_obstacle_potential|rpo_chomp_obstacle_potential]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:176-176`
<!-- vulcan:connections:end -->

## Limitations
Each iteration allocates a column slice; the clearance query dominates runtime and is called once per sample with no caching between finite-difference evaluations, which makes `rpo_chomp_numeric_gradient` expensive for large `n_waypoints`. Averaging means a short deep penetration on a long path is diluted relative to a sum-based penalty.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 169.
