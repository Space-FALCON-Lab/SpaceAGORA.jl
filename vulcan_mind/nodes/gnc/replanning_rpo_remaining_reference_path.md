---
id: gnc.replanning_rpo_remaining_reference_path
label: rpo_remaining_reference_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_remaining_reference_path
  lines:
  - 172
  - 172
inputs:
- id: plan
  type: RPOPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: current_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `current_rtn`.
- id: sample_ds_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `sample_ds_m` (default `0.10`).
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
  description: Return value of `rpo_remaining_reference_path`. Returns `reshape(collect(current),
    3, 1)` or `rpo_sample_path_polyline(path, sample_ds_m)`.
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

# rpo_remaining_reference_path

## Purpose
Extracts the portion of the current plan's reference trajectory that still lies ahead of the chaser, starting from the chaser's actual position, and resamples it at a uniform arc-length spacing so clearance statistics and retiming operate on the path yet to be flown rather than the whole plan.

## Design & Implementation
Signature `rpo_remaining_reference_path(plan::RPOPlan, current_rtn; sample_ds_m::Real=0.10)`. If `plan.r_ref_rtn` has zero columns it returns the current position as a `3 x 1` matrix. Otherwise it scans every column with an `@inbounds` loop to find the index minimising squared Euclidean distance to `current`, takes `tail = plan.r_ref_rtn[:, best_idx:end]`, prepends `current` as the first column of a new `3 x (length(tail)+1)` matrix, and returns `rpo_sample_path_polyline(path, sample_ds_m)`, the polyline resampled at `sample_ds_m` metre intervals.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RPOPlan | n/a | yes | Positional argument `plan`. |
| in | `current_rtn` | Any | n/a | yes | Positional argument `current_rtn`. |
| in | `sample_ds_m` | Real | n/a | no | Keyword argument `sample_ds_m` (default `0.10`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_remaining_reference_path`. Returns `reshape(collect(current), 3, 1)` or `rpo_sample_path_polyline(path, sample_ds_m)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_retime_existing_plan|rpo_retime_existing_plan]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:265-265`
- [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:238-238`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:193-193`
<!-- vulcan:connections:end -->

## Limitations
Nearest-point selection is purely spatial, so on a path that loops back near itself the tail may jump to a later or earlier segment than the one actually being tracked; no time or progress information from `plan.t_ref_s` is used. When the nearest reference point is behind the chaser, the segment from `current` back to that point is included and doubles back. The linear scan is O(n) per call with no caching of the previous index.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 172.
