---
id: analysis.telemetry_loading__initial_time_et
label: _initial_time_et
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _initial_time_et
  lines:
  - 22
  - 22
inputs:
- id: initial_time
  type: InitialTime
  units: n/a
  required: true
  description: Positional argument `initial_time`.
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
  description: Return value of `_initial_time_et`.
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

# _initial_time_et

## Purpose
Converts a scenario `InitialTime` (calendar fields) into SPICE ephemeris time (TDB seconds past J2000) so telemetry samples offset in seconds from the scenario start can be transformed between inertial and planet-fixed frames at the correct epoch.

## Design & Implementation
Formats the fields with `@sprintf("%04d-%02d-%02dT%02d:%02d:%09.6f", year, month, day, hour, minute, second)`, converting integer fields with `Int` and the seconds with `Float64` so fractional seconds are preserved to microseconds, then calls `utc2et` on the ISO-style string. The return type annotation forces `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | InitialTime | n/a | yes | Positional argument `initial_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_initial_time_et`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:148-148`
- [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:105-105`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
`utc2et` requires a leap-seconds kernel (`naif*.tls`) to be furnished; otherwise it throws a SPICE error that is not caught here. Seconds of `60` or greater format as a legal-looking string that SPICE rejects. Years before 1000 or after 9999 break the `%04d` layout. The function is called per scenario, not per sample, so cost is negligible, but it assumes `InitialTime` fields are UTC, not TDB or TT.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 22.
