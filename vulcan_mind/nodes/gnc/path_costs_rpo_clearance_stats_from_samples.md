---
id: gnc.path_costs_rpo_clearance_stats_from_samples
label: rpo_clearance_stats_from_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_clearance_stats_from_samples
  lines:
  - 31
  - 31
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
  type: Real
  units: n/a
  required: true
  description: Positional argument `safe_distance_m`.
- id: cost_cutoff
  type: Real
  units: n/a
  required: false
  description: Keyword argument `cost_cutoff` (default `Inf`).
- id: w_obs
  type: Real
  units: n/a
  required: false
  description: Keyword argument `w_obs` (default `0.0`).
- id: obstacle_sigmoid_k
  type: Real
  units: n/a
  required: false
  description: Keyword argument `obstacle_sigmoid_k` (default `1.0e6`).
- id: obstacle_sigmoid_tol_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `obstacle_sigmoid_tol_m` (default `0.0`).
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
  description: Return value of `rpo_clearance_stats_from_samples`. Returns `(`.
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

# rpo_clearance_stats_from_samples

## Purpose
`rpo_clearance_stats_from_samples` is the collision-assessment sweep of the RPO cost function. It walks every sampled point of a candidate path, queries its clearance against the station geometry, and returns the worst-case clearance, the count and fraction of samples in violation, and the accumulated smooth obstacle score that feeds the weighted objective — with an early-exit so hopeless candidates are abandoned partway through.

## Design & Implementation
State is three accumulators: `min_clearance = Inf`, `violation_count = 0` and `obstacle_score = 0.0`. `threshold` is precomputed once as `safe_distance_m - obstacle_sigmoid_tol_m`, `k` as `Float64(obstacle_sigmoid_k)`. The `@inbounds` loop over `j in 1:size(samples, 2)` builds an `SVector{3, Float64}` from the column, calls `rpo_clearance_distance_to_station(p, geometry)`, folds the result into `min_clearance`, adds `rpo_obstacle_sigmoid_penalty(clearance, threshold, k)` to the score, and increments `violation_count` when `clearance < 0.0` — note that violations are counted against zero clearance, not against `safe_distance_m`. After each sample, if `w_obs > 0.0`, `cost_cutoff` is finite and `w_obs * obstacle_score > cost_cutoff`, it returns immediately with `cutoff_exceeded = true`. Both exits return the same five-field NamedTuple, with `violation_fraction = violation_count / max(size(samples, 2), 1)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `cost_cutoff` | Real | n/a | no | Keyword argument `cost_cutoff` (default `Inf`). |
| in | `w_obs` | Real | n/a | no | Keyword argument `w_obs` (default `0.0`). |
| in | `obstacle_sigmoid_k` | Real | n/a | no | Keyword argument `obstacle_sigmoid_k` (default `1.0e6`). |
| in | `obstacle_sigmoid_tol_m` | Real | n/a | no | Keyword argument `obstacle_sigmoid_tol_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_clearance_stats_from_samples`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:101-101`
- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:60-60`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:43-43`
- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:48-48`
- `callees` → [[gnc.path_costs_rpo_obstacle_sigmoid_penalty|rpo_obstacle_sigmoid_penalty]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:50-50`
<!-- vulcan:connections:end -->

## Limitations
On the early-exit path `violation_fraction` is divided by the *total* sample count even though only a prefix was examined, so the reported fraction understates the true violation rate — and `min_clearance` likewise reflects only the prefix. A path evaluated with zero samples returns `min_clearance = Inf` and zero violations, i.e. it is reported as perfectly safe. Because the sigmoid saturates at the default sharpness, `obstacle_score` is effectively an integer count of sub-threshold samples rather than a graded measure. Clearance is sampled pointwise, so an obstacle thinner than the sampling step between two consecutive points is not detected at all — tunnelling is possible. The cutoff check runs every iteration, adding a branch to the innermost loop of the planner.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_costs.jl` line 31.
