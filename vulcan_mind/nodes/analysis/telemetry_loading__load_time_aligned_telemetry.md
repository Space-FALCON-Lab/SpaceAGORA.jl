---
id: analysis.telemetry_loading__load_time_aligned_telemetry
label: _load_time_aligned_telemetry
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _load_time_aligned_telemetry
  lines:
  - 213
  - 213
inputs:
- id: cfg
  type: TimeAlignedScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: max_points
  type: Int
  units: n/a
  required: true
  description: Positional argument `max_points`.
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
  description: Return value of `_load_time_aligned_telemetry`. Returns `(`.
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

# _load_time_aligned_telemetry

## Purpose
Loads and normalises a time-aligned telemetry Arrow file for one `TimeAlignedScenarioConfig`: it extracts time, altitude and J2000 position (km), optional per-sample velocity truth, optional Keplerian or Cartesian initial-condition columns, sorts by time, truncates, rebases time to zero, applies an optional day/night illumination mask and returns a flat named tuple consumed by the verification runners.

## Design & Implementation
Every column name comes from `cfg` and is fetched through `_require_column`. The six Keplerian IC columns (`sma`, `ecc`, `inc`, `aop`, `raan`, `ta`) and six Cartesian IC columns must be all-or-nothing, and at least one full set must exist; the three velocity-truth columns are likewise all-or-nothing. Violations throw `ArgumentError` naming `cfg.name`. Rows are reordered by `sortperm(time_s)`; `max_points > 0` keeps the first `max_points`; at least two samples are required and `time_s .- t0` rebases. Keplerian IC values are taken from the first sorted row before masking; Cartesian IC values from `perm[1]` of the unsorted columns. If `cfg.truth_mask` is `:dayside` or `:nightside`, `comparison_frame` must be `:inertial`; each sample is classified by the sign of `dot(r, s_hat)` with `s_hat` from `_sun_unit_vector_j2000`, a summary line is printed, and fewer than two survivors throws. Velocities are the truth columns when present, otherwise `_differentiate_series` of each position component.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | TimeAlignedScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `max_points` | Int | n/a | yes | Positional argument `max_points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_load_time_aligned_telemetry`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:219-219`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:382-382`
- `callees` → [[analysis.telemetry_loading__differentiate_series|_differentiate_series]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:373-373`
- `callees` → [[analysis.telemetry_loading__require_column|_require_column]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:216-216`
- `callees` → [[analysis.telemetry_loading__sun_unit_vector_j2000|_sun_unit_vector_j2000]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:338-338`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:345-345`
<!-- vulcan:connections:end -->

## Limitations
The illumination screen is a pure geometric half-space test (`cosang > 0`), ignoring the planet's shadow cylinder, so samples in eclipse behind the planet are classed as lit. Absent IC columns are represented as `NaN` vectors sized `n_rows` and then sorted and sliced needlessly. The mask applies after truncation, so `max_points` and `truth_mask` interact and can leave very few samples. Differentiated velocities from position telemetry amplify noise and are first-order at the ends. The function prints to stdout unconditionally when a mask is active.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 213.
