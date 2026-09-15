---
id: gnc.pso_adaptive_policy_rpo_estimate_geometry_complexity
label: rpo_estimate_geometry_complexity
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl
  symbol: rpo_estimate_geometry_complexity
  lines:
  - 2
  - 2
inputs:
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
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: sample_ds_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `sample_ds_m` (default `0.25`).
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
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
  description: Return value of `rpo_estimate_geometry_complexity`. Returns `clamp(0.7
    * buffer_fraction + 0.3 * clearance_term, 0.0, 1.0)`.
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

# rpo_estimate_geometry_complexity

## Purpose
Scores how hard a straight-line approach between two RTN points would be, so the planner can decide whether to spend effort on particle-swarm path search.

## Theory & Math
$$
c = \mathrm{clamp}\left( 0.7 \, b + 0.3 \, \frac{1}{1 + d_{\min}},\; 0,\; 1 \right)
$$

with $d_{\min}$ the minimum clearance in metres and $b \in \{0,1\}$ indicating whether $d_{\min}$ falls below the safe distance.

## Design & Implementation
Samples the straight segment from start to goal at `sample_ds_m` spacing and measures clearance statistics along it. Two terms combine: `buffer_fraction` is one when the minimum clearance falls below the safe distance, and `clearance_term` is one when clearance is non-positive, otherwise the reciprocal `1/(1+d)` so complexity decays smoothly as clearance grows. The weighted sum, 0.7 on the buffer term and 0.3 on clearance, is clamped to the unit interval.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `sample_ds_m` | Real | n/a | no | Keyword argument `sample_ds_m` (default `0.25`). |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_estimate_geometry_complexity`. Returns `clamp(0.7 * buffer_fraction + 0.3 * clearance_term, 0.0, 1.0)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_adaptive_policy_rpo_adaptive_pso_config|rpo_adaptive_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:31-31`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:7-7`
- `callees` → [[gnc.clearance_rpo_path_clearance_stats|rpo_path_clearance_stats]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:5-5`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:4-4`
<!-- vulcan:connections:end -->

## Limitations
Only the straight line is probed, so a case where the direct path is clear but the reachable corridor is not scores as easy; the 0.7 and 0.3 weights are fixed constants with no configuration hook.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl` line 2.
