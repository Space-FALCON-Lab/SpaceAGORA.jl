---
id: gnc.pso_adaptive_policy_rpo_probe_geometry_metrics
label: rpo_probe_geometry_metrics
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl
  symbol: rpo_probe_geometry_metrics
  lines:
  - 12
  - 12
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
  description: Return value of `rpo_probe_geometry_metrics`. Returns `(`.
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

# rpo_probe_geometry_metrics

## Purpose
Measures the straight-line approach between two RTN points and reports whether it is flyable, forming the baseline that searched paths are compared against.

## Design & Implementation
Builds the two-point segment, samples it at `sample_ds_m`, and evaluates clearance statistics. It returns a named tuple carrying `min_clearance` and `violation_fraction` straight through, a `detour_ratio` fixed at 1.0 because the straight path is by definition its own reference, and a `success` flag testing minimum clearance against the safe distance with a 1e-9 tolerance to absorb floating-point equality at the boundary.

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
| out | `result` | Any | n/a | — | Return value of `rpo_probe_geometry_metrics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:20-20`
- `callees` → [[gnc.clearance_rpo_path_clearance_stats|rpo_path_clearance_stats]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:15-15`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:14-14`
<!-- vulcan:connections:end -->

## Limitations
The hard-coded unit detour ratio means this probe cannot express path lengthening; it is only meaningful as the denominator other candidate paths are measured against.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl` line 12.
