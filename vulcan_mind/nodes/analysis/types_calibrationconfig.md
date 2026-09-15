---
id: analysis.types_calibrationconfig
label: CalibrationConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: CalibrationConfig
  lines:
  - 44
  - 44
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `false`).
- id: profiles
  type: Vector{Symbol}
  units: n/a
  required: false
  description: Field `profiles` (default `Symbol[:full]`).
- id: search_on_quick_subset
  type: Bool
  units: n/a
  required: false
  description: Field `search_on_quick_subset` (default `true`).
- id: fit_cd_scale
  type: Bool
  units: n/a
  required: false
  description: Field `fit_cd_scale` (default `true`).
- id: cd_scale_min
  type: Float64
  units: n/a
  required: false
  description: Field `cd_scale_min` (default `0.85`).
- id: cd_scale_max
  type: Float64
  units: n/a
  required: false
  description: Field `cd_scale_max` (default `1.15`).
- id: cd_scale_steps
  type: Int
  units: n/a
  required: false
  description: Field `cd_scale_steps` (default `3`).
- id: fit_cr
  type: Bool
  units: n/a
  required: false
  description: Field `fit_cr` (default `true`).
- id: cr_min
  type: Float64
  units: n/a
  required: false
  description: Field `cr_min` (default `1.15`).
- id: cr_max
  type: Float64
  units: n/a
  required: false
  description: Field `cr_max` (default `1.45`).
- id: cr_steps
  type: Int
  units: n/a
  required: false
  description: Field `cr_steps` (default `3`).
- id: fit_bias
  type: Bool
  units: n/a
  required: false
  description: Field `fit_bias` (default `true`).
- id: bias_abs_max_km
  type: Float64
  units: n/a
  required: false
  description: Field `bias_abs_max_km` (default `500.0`).
- id: objective
  type: String
  units: n/a
  required: false
  description: Field `objective` (default `"mean_nmae"`).
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
  type: CalibrationConfig
  units: n/a
  description: Constructed `CalibrationConfig` (keyword constructor via @kwdef).
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

# CalibrationConfig

## Purpose
Keyword struct controlling the optional parameter-calibration sweep in telemetry verification, where drag-coefficient scale, SRP reflectivity `Cr`, and an along-track bias are grid-searched to minimise an error objective before scoring.

## Design & Implementation
`enabled=false` gates the whole sweep; `profiles=[:full]` lists which study profiles calibrate and `search_on_quick_subset=true` lets the grid search run on the smaller quick comparison set. Each fitted quantity has an enable flag and range: `fit_cd_scale` over `[cd_scale_min=0.85, cd_scale_max=1.15]` in `cd_scale_steps=3` points, `fit_cr` over `[cr_min=1.15, cr_max=1.45]` in `cr_steps=3` points, and `fit_bias` bounded by `bias_abs_max_km=500.0`. `objective="mean_nmae"` names the scalar metric minimised. All fields are immutable after construction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `false`). |
| in | `profiles` | Vector{Symbol} | n/a | no | Field `profiles` (default `Symbol[:full]`). |
| in | `search_on_quick_subset` | Bool | n/a | no | Field `search_on_quick_subset` (default `true`). |
| in | `fit_cd_scale` | Bool | n/a | no | Field `fit_cd_scale` (default `true`). |
| in | `cd_scale_min` | Float64 | n/a | no | Field `cd_scale_min` (default `0.85`). |
| in | `cd_scale_max` | Float64 | n/a | no | Field `cd_scale_max` (default `1.15`). |
| in | `cd_scale_steps` | Int | n/a | no | Field `cd_scale_steps` (default `3`). |
| in | `fit_cr` | Bool | n/a | no | Field `fit_cr` (default `true`). |
| in | `cr_min` | Float64 | n/a | no | Field `cr_min` (default `1.15`). |
| in | `cr_max` | Float64 | n/a | no | Field `cr_max` (default `1.45`). |
| in | `cr_steps` | Int | n/a | no | Field `cr_steps` (default `3`). |
| in | `fit_bias` | Bool | n/a | no | Field `fit_bias` (default `true`). |
| in | `bias_abs_max_km` | Float64 | n/a | no | Field `bias_abs_max_km` (default `500.0`). |
| in | `objective` | String | n/a | no | Field `objective` (default `"mean_nmae"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CalibrationConfig | n/a | — | Constructed `CalibrationConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:390-390`
- [[analysis.types_orbiteventsscenarioconfig|OrbitEventsScenarioConfig]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:152-152`
- [[analysis.types_timealignedscenarioconfig|TimeAlignedScenarioConfig]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:215-215`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Step counts of 3 give a very coarse grid; there is no refinement stage described by the struct. Ranges are not validated (`cd_scale_min > cd_scale_max` or `cd_scale_steps < 1` are accepted). The objective is a free-form string, so an unsupported name fails only when the calibration routine looks it up. `bias_abs_max_km` is in kilometres while most other scenario fields use metres, which is an easy unit slip for manifest authors.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 44.
