---
id: analysis.telemetry_loading__load_telemetry_curve
label: _load_telemetry_curve
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _load_telemetry_curve
  lines:
  - 200
  - 200
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
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
  description: Return value of `_load_telemetry_curve`. Returns `(`.
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

# _load_telemetry_curve

## Purpose
Loads a pre-reduced telemetry curve (one altitude value per orbit number) from an Arrow file for scenarios that compare periapsis or apoapsis altitude by orbit index rather than by time.

## Design & Implementation
Reads `DataFrame(Arrow.Table(path))`, asserts that both `orbit` and `altitude` columns exist, and sorts in place by `:orbit`. When `max_points > 0` and the frame is longer, it keeps only `first(df, max_points)` rows. Returns a named tuple `(orbit, altitude)` of `Vector{Float64}` built through `_to_float_vector` with contexts `"telemetry-orbit"` and `"telemetry-altitude"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `max_points` | Int | n/a | yes | Positional argument `max_points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_load_telemetry_curve`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_error_tables_orbit_rows_errors|_orbit_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:9-9`

**Downstream**

- `callees` → [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:208-208`
<!-- vulcan:connections:end -->

## Limitations
The column presence check uses `@assert`, which is compiled out under `--check-bounds=no`/`-O3` with assertions disabled and then produces a confusing `ArgumentError` from `getproperty` instead. Units of `altitude` are not validated (the rest of the pipeline assumes km). Truncation keeps the earliest orbits only, with no option to subsample. Duplicate orbit numbers survive the sort and are passed through.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 200.
