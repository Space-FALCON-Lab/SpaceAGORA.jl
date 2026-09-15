---
id: analysis.calibration__annotate_calibration_rows
label: _annotate_calibration_rows
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _annotate_calibration_rows
  lines:
  - 70
  - 70
inputs:
- id: rows
  type: AbstractVector{<:NamedTuple}
  units: n/a
  required: true
  description: Positional argument `rows`.
- id: cd_scale
  type: Float64
  units: n/a
  required: true
  description: Positional argument `cd_scale`.
- id: cr_value
  type: Float64
  units: n/a
  required: true
  description: Positional argument `cr_value`.
- id: bias_by_event
  type: Dict{String, Float64}
  units: n/a
  required: true
  description: Positional argument `bias_by_event`.
- id: score
  type: Float64
  units: n/a
  required: true
  description: Positional argument `score`.
- id: selected_runtime_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `selected_runtime_s`.
- id: dt_max_orbit_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_max_orbit_s`.
- id: calibration_runtime_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `calibration_runtime_s`.
- id: calibration_used
  type: Bool
  units: n/a
  required: true
  description: Positional argument `calibration_used`.
- id: solver_info
  type: Any
  units: n/a
  required: true
  description: Positional argument `solver_info`.
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
  type: Vector{NamedTuple}
  units: n/a
  description: Return value of `_annotate_calibration_rows`.
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

# _annotate_calibration_rows

## Purpose

`_annotate_calibration_rows` stamps the calibration and solver provenance onto every output row of a telemetry-verification run, so the emitted table records not just the errors but exactly which drag scale, reflectivity, per-event bias and solver path produced them.

## Design & Implementation

The function iterates the input `rows` and `merge`s each `NamedTuple` with a fixed block of eighteen additional fields, pushing the results into a `Vector{NamedTuple}`. Scalars `calibration_used`, `calibrated_cd_scale`, `calibrated_cr`, `calibration_score`, `selected_simulation_runtime_s`, `dt_max_orbit_s` and `calibration_runtime_s` are copied verbatim onto every row, while `calibrated_bias_km` is looked up per row with `get(bias_by_event, String(row.event), 0.0)` so an event with no fitted bias records a clean zero. Eight further fields are read off `solver_info`, covering solver mode, sequence, fallback usage and count, fallback trigger, return code, `maxiters` and whether a `maxiters` retry was used. The input rows are not mutated; a new vector is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rows` | AbstractVector{<:NamedTuple} | n/a | yes | Positional argument `rows`. |
| in | `cd_scale` | Float64 | n/a | yes | Positional argument `cd_scale`. |
| in | `cr_value` | Float64 | n/a | yes | Positional argument `cr_value`. |
| in | `bias_by_event` | Dict{String, Float64} | n/a | yes | Positional argument `bias_by_event`. |
| in | `score` | Float64 | n/a | yes | Positional argument `score`. |
| in | `selected_runtime_s` | Float64 | n/a | yes | Positional argument `selected_runtime_s`. |
| in | `dt_max_orbit_s` | Float64 | n/a | yes | Positional argument `dt_max_orbit_s`. |
| in | `calibration_runtime_s` | Float64 | n/a | yes | Positional argument `calibration_runtime_s`. |
| in | `calibration_used` | Bool | n/a | yes | Positional argument `calibration_used`. |
| in | `solver_info` | Any | n/a | yes | Positional argument `solver_info`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{NamedTuple} | n/a | — | Return value of `_annotate_calibration_rows`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:197-197`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations

The returned container is `Vector{NamedTuple}` without a concrete element type, so downstream code loses type stability and pays dynamic dispatch when consuming it. The eighteen appended field names are hard-coded, meaning any change to the solver telemetry schema requires editing this function and every reader. `merge` lets an appended field silently overwrite an identically named field already on the row, and `solver_info` is accessed by property without any check that all eight properties exist.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/calibration.jl` line 70.
