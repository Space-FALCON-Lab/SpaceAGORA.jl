---
id: analysis.calibration__grid_values
label: _grid_values
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _grid_values
  lines:
  - 23
  - 23
inputs:
- id: min_v
  type: Float64
  units: n/a
  required: true
  description: Positional argument `min_v`.
- id: max_v
  type: Float64
  units: n/a
  required: true
  description: Positional argument `max_v`.
- id: steps
  type: Int
  units: n/a
  required: true
  description: Positional argument `steps`.
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
  type: Vector{Float64}
  units: n/a
  description: Return value of `_grid_values`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _grid_values

## Purpose

`_grid_values(min_v, max_v, steps)` builds the candidate vector swept by the calibration search, returning a `Vector{Float64}` of `steps` values spanning `min_v` to `max_v` inclusive. It is used for both the drag-scale and the reflectivity coefficient grids.

## Design & Implementation

The function collapses two degenerate cases to a single-element grid `[min_v]`: when `steps <= 1`, and when the endpoints are equal within `isapprox(min_v, max_v; rtol=0.0, atol=1e-12)`. Note that the tolerance is purely absolute, `rtol` being forced to zero, so the collapse triggers only for endpoints within 1e-12 of each other regardless of magnitude. Otherwise it returns `collect(range(min_v, max_v, length=steps))`, a uniformly spaced inclusive grid materialised eagerly as a dense vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `min_v` | Float64 | n/a | yes | Positional argument `min_v`. |
| in | `max_v` | Float64 | n/a | yes | Positional argument `max_v`. |
| in | `steps` | Int | n/a | yes | Positional argument `steps`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Float64} | n/a | — | Return value of `_grid_values`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:149-149`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

No check enforces `min_v <= max_v`; a reversed pair yields a descending grid rather than an error. The absolute 1e-12 tolerance is inappropriate for parameters whose natural scale is far from unity, and the spacing is always linear, so a parameter better searched logarithmically must be transformed by the caller. Large `steps` values allocate the entire grid up front.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/calibration.jl` line 23.
