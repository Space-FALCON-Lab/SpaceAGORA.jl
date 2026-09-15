---
id: gnc.clearance_rpo_path_clearance_stats
label: rpo_path_clearance_stats
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/clearance.jl
  symbol: rpo_path_clearance_stats
  lines:
  - 15
  - 15
inputs:
- id: path_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `path_body`.
- id: geometry
  type: RPOReferenceGeometry
  units: n/a
  required: true
  description: Positional argument `geometry`.
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
  description: Return value of `rpo_path_clearance_stats`. Returns `(`.
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

# rpo_path_clearance_stats

## Purpose
Summarises how safely a candidate RPO path passes the station, giving planners a single set of numbers to compare trajectories.

## Design & Implementation
Validates that `path_body` is a 3 by N matrix and raises `ArgumentError` otherwise. It then walks every column, evaluating `rpo_clearance_distance_to_station` at each sample, tracking the running minimum and counting samples whose clearance falls below `safe_distance_m`. Returns a named tuple of `min_clearance`, `violation_count` and `violation_fraction`, the last guarded by `max(N, 1)` so an empty path cannot divide by zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path_body` | Any | n/a | yes | Positional argument `path_body`. |
| in | `geometry` | RPOReferenceGeometry | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_path_clearance_stats`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_adaptive_policy_rpo_estimate_geometry_complexity|rpo_estimate_geometry_complexity]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:5-5`
- [[gnc.pso_adaptive_policy_rpo_probe_geometry_metrics|rpo_probe_geometry_metrics]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:15-15`
- [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:239-239`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:18-18`
- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
Clearance is evaluated only at supplied samples, so a path that dips inside the keep-out between two widely spaced waypoints is scored as clear; sampling density is the caller's responsibility.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/clearance.jl` line 15.
