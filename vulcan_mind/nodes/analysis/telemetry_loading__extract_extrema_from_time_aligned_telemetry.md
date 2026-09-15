---
id: analysis.telemetry_loading__extract_extrema_from_time_aligned_telemetry
label: _extract_extrema_from_time_aligned_telemetry
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _extract_extrema_from_time_aligned_telemetry
  lines:
  - 429
  - 429
inputs:
- id: time_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `time_s`.
- id: altitude_km
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `altitude_km`.
- id: min_sep_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `min_sep_s`.
- id: speed_kmps
  type: Union{Nothing, Vector{Float64}}
  units: n/a
  required: false
  description: Keyword argument `speed_kmps` (default `nothing`).
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
  description: Return value of `_extract_extrema_from_time_aligned_telemetry`. Returns
    `(`.
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

# _extract_extrema_from_time_aligned_telemetry

## Purpose
Finds periapsis and apoapsis events in a telemetry altitude history by locating strict local minima and maxima, de-duplicated by a minimum separation, so that time-aligned telemetry can be compared with simulation extrema on a per-orbit basis.

## Design & Implementation
Arguments are `time_s`, `altitude_km` (equal length, `n >= 3`) and `min_sep_s`; an optional `speed_kmps` vector of the same length provides the speed at each event, else `NaN`. The `@inbounds` loop over `2:(n-1)` classifies sample `i` as a periapsis when `a1 <= a0 && a1 < a2` and an apoapsis when `a1 >= a0 && a1 > a2`. A new event is pushed if the list is empty or `ti - last_t >= min_sep_s`; otherwise the last event is replaced only when the new one is more extreme. Empty periapsis or apoapsis lists throw `ArgumentError`. Returns `(peri=(orbit, altitude, time_s, speed_kmps), apo=(...))` with `orbit = collect(1.0:count)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `time_s` | Vector{Float64} | n/a | yes | Positional argument `time_s`. |
| in | `altitude_km` | Vector{Float64} | n/a | yes | Positional argument `altitude_km`. |
| in | `min_sep_s` | Float64 | n/a | yes | Positional argument `min_sep_s`. |
| in | `speed_kmps` | Union{Nothing, Vector{Float64}} | n/a | no | Keyword argument `speed_kmps` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_extract_extrema_from_time_aligned_telemetry`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:85-85`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:464-464`
<!-- vulcan:connections:end -->

## Limitations
Events are taken at sample times without interpolation, so the reported extremum is biased upward for periapsis and downward for apoapsis by up to the altitude change over half a sample interval. The asymmetric `<=`/`<` test registers the last point of a flat run, and a strictly monotone history yields no events and throws. `min_sep_s` is applied separately to periapsis and apoapsis, so a periapsis and apoapsis may be arbitrarily close. Orbit numbering starts at 1 regardless of where the telemetry begins, so it may be offset from the simulation's orbit counter.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 429.
