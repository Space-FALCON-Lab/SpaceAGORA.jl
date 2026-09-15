---
id: analysis.calibration__calibration_score
label: _calibration_score
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _calibration_score
  lines:
  - 58
  - 58
inputs:
- id: rows
  type: AbstractVector{<:NamedTuple}
  units: n/a
  required: true
  description: Positional argument `rows`.
- id: objective
  type: String
  units: n/a
  required: true
  description: Positional argument `objective`.
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
  type: Float64
  units: n/a
  description: Return value of `_calibration_score`.
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

# _calibration_score

## Purpose

`_calibration_score(rows, objective)` reduces a set of per-event error rows to the single scalar the calibration grid search minimises. It supports three objectives selected by the `objective` string: `mean_nmae`, `mean_rmse_km` and `max_nmae`.

## Design & Implementation

An empty `rows` vector short-circuits to `Inf`, so a candidate that produced no comparable events can never win the search. Otherwise the function extracts the relevant field from every `NamedTuple` row into a `Vector{Float64}` and reduces it: `mean` over `r.nmae` for `mean_nmae`, `mean` over `r.rmse_km` for `mean_rmse_km`, and `maximum` over `r.nmae` for `max_nmae`. Any other string falls through to `throw(ArgumentError("Unsupported calibration objective '$objective'"))`, so a typo in configuration fails loudly rather than silently defaulting.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rows` | AbstractVector{<:NamedTuple} | n/a | yes | Positional argument `rows`. |
| in | `objective` | String | n/a | yes | Positional argument `objective`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_calibration_score`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:171-171`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:61-61`
<!-- vulcan:connections:end -->

## Limitations

Objective selection is by string comparison at every call rather than by dispatch, and each branch builds an intermediate array before reducing, so the function allocates once per candidate evaluation. Rows are weighted equally: a short event counts as much as a long one, and no per-event weighting or outlier rejection is available. `NaN` in any row propagates into the mean, quietly poisoning the comparison since `NaN` loses every ordering test.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/calibration.jl` line 58.
