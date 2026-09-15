---
id: analysis.ic_fit__ic_fit_mean_axis_rmse
label: _ic_fit_mean_axis_rmse
kind: function
source:
  file: src/analysis/verification/telemetry_verification/ic_fit.jl
  symbol: _ic_fit_mean_axis_rmse
  lines:
  - 66
  - 66
inputs:
- id: summary
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `summary`.
- id: scenario_name
  type: String
  units: n/a
  required: true
  description: Positional argument `scenario_name`.
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
  description: Return value of `_ic_fit_mean_axis_rmse`.
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

# _ic_fit_mean_axis_rmse

## Purpose

Collapses a verification summary table into the single scalar the fit reports as its before/after quality metric: the mean of the x, y and z position RMSE values, in kilometres, for one scenario.

## Design & Implementation

Filters `summary::DataFrame` to rows matching `scenario_name` whose `event` lies in `_IC_FIT_EVENTS`, insists on exactly three surviving rows (otherwise `ArgumentError`, catching a run that silently lost an axis), and returns `sum(Float64.(rows.rmse_km)) / 3.0`. The `nrow(rows) == 3` check is the guard that keeps a partially failed propagation from producing an optimistic-looking metric.

## Theory & Math

With per-axis root-mean-square errors $\varepsilon_x, \varepsilon_y, \varepsilon_z$ over the comparison samples, the reported metric is the arithmetic mean

$$\bar{\varepsilon} = \frac{\varepsilon_x + \varepsilon_y + \varepsilon_z}{3}$$

which is not the RMS norm of the position error vector $\sqrt{\varepsilon_x^2 + \varepsilon_y^2 + \varepsilon_z^2}$; it is smaller by construction and only proportional to it when the three axes carry equal error.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `summary` | DataFrame | n/a | yes | Positional argument `summary`. |
| in | `scenario_name` | String | n/a | yes | Positional argument `scenario_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_ic_fit_mean_axis_rmse`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_ic_fit_fit_initial_state|fit_initial_state]] · `callees` → `callers` · feedback · `src/analysis/verification/telemetry_verification/ic_fit.jl:165-165`

**Downstream**

- `callees` → [[envana.ana_ic_fit_fit_initial_state|fit_initial_state]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations

The equal-weight average hides axis anisotropy, so an along-track blow-up on one axis can be masked by two well-fitted axes. The hard requirement of exactly three rows means an extra duplicated event row aborts the fit rather than being deduplicated. The `rmse_km` column is assumed present and finite; a `NaN` from a diverged run propagates straight into the returned mean without being flagged.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/ic_fit.jl` line 66.
