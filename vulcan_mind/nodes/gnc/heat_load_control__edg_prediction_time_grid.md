---
id: gnc.heat_load_control__edg_prediction_time_grid
label: _edg_prediction_time_grid
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_prediction_time_grid
  lines:
  - 105
  - 105
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
  description: Return value of `_edg_prediction_time_grid`. Returns `collect(range(0.0,
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

# _edg_prediction_time_grid

## Purpose
Builds the uniformly spaced time grid over which heat-load trajectories, costates, and switching profiles are evaluated.

## Design & Implementation
Takes `duration::Float64` in seconds and a fixed nominal `step = 1.0` s. The point count is `clamp(ceil(Int, duration / step) + 1, 64, 2000)`, and the result is `collect(range(0.0, duration; length=n))`, so the actual spacing is `duration / (n - 1)` and may differ from 1 s at either clamp bound. Returns a `Vector{Float64}` beginning at 0.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `duration` | Float64 | n/a | yes | Positional argument `duration`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_prediction_time_grid`. Returns `collect(range(0.0, duration; length=n))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:135-135`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Passes longer than 1999 s are sampled coarser than 1 s, degrading the accuracy of the explicit Euler costate integration and the RK4 trajectory. Very short passes are oversampled to 64 points. The grid is always uniform, with no refinement near periapsis where heating peaks.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 105.
