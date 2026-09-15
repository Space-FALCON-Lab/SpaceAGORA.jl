---
id: gnc.replanning_rpo_reference_tracking_error
label: rpo_reference_tracking_error
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_reference_tracking_error
  lines:
  - 197
  - 197
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
  description: Return value of `rpo_reference_tracking_error`. Returns `best_dist`.
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

# rpo_reference_tracking_error

## Purpose
Measures how far the chaser has drifted from the planned reference trajectory, returning the minimum Euclidean distance from the current RTN position to any reference sample. `rpo_replanning_decision` compares this against `tracking_error_retime_m` and `tracking_error_replan_m` when no obstacles are active.

## Design & Implementation
Signature `rpo_reference_tracking_error(plan::RPOPlan, current_rtn)`. It converts `current_rtn` to an `SVector{3,Float64}`, returns `Inf` if `plan.r_ref_rtn` has no columns, and otherwise iterates all columns with `@inbounds`, computing `sqrt(dx^2 + dy^2 + dz^2)` for each and keeping the minimum with `min`. The result is in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RPOPlan | n/a | yes | Positional argument `plan`. |
| in | `current_rtn` | Any | n/a | yes | Positional argument `current_rtn`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_reference_tracking_error`. Returns `best_dist`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:229-229`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Distance is to the discrete sample points, not to the polyline between them, so with coarse reference spacing the reported error is biased upward by up to half the sample spacing. `sqrt` is evaluated for every column rather than comparing squared distances once, a minor cost in a hot loop. The metric ignores time, so a chaser that is exactly on the path but far behind schedule registers zero error.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 197.
