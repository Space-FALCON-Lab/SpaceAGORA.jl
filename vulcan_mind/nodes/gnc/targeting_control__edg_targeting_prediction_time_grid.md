---
id: gnc.targeting_control__edg_targeting_prediction_time_grid
label: _edg_targeting_prediction_time_grid
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_prediction_time_grid
  lines:
  - 409
  - 409
inputs:
- id: duration
  type: Float64
  units: n/a
  required: true
  description: Positional argument `duration`.
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
  description: Return value of `_edg_targeting_prediction_time_grid`. Returns `collect(range(0.0,
    duration; length=n))`.
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

# _edg_targeting_prediction_time_grid

## Purpose
Builds the uniform time grid the prediction integrates over for one drag passage.

## Design & Implementation
Targets a 0.1 s step, computes the point count as the ceiling of `duration / 0.1` plus one, clamps it between 64 and 20,000, and returns `range(0, duration; length=n)` collected.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `duration` | Float64 | n/a | yes | Positional argument `duration`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_targeting_prediction_time_grid`. Returns `collect(range(0.0, duration; length=n))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:726-726`
- [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:684-684`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For passages longer than 2,000 s the clamp coarsens the step above 0.1 s without warning; for very short passages the 64-point floor makes the step finer than needed.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 409.
