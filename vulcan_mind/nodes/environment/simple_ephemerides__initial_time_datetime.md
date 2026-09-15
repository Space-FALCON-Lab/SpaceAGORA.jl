---
id: environment.simple_ephemerides__initial_time_datetime
label: _initial_time_datetime
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _initial_time_datetime
  lines:
  - 40
  - 40
inputs:
- id: initial_time
  type: Any
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
  type: DateTime
  units: n/a
  description: Return value of `_initial_time_datetime`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _initial_time_datetime

## Purpose
Converts a mission `initial_time` record (with separate `year`, `month`, `day`, `hour`, `minute`, and fractional `second` fields) into a `Dates.DateTime`, the common starting point for both SPICE and simple ephemeris time conversions.

## Design & Implementation
A `DateTime` is constructed from the integer parts `Int(initial_time.year)` through `Int(initial_time.minute)` with the seconds field set to 0. The fractional seconds are then added as `Millisecond(round(Int, 1000 * Float64(initial_time.second)))`, so sub-second precision is retained to the nearest millisecond. The function is `@inline` and duck-typed on `initial_time`, accepting any object exposing the six fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DateTime | n/a | — | Return value of `_initial_time_datetime`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:53-53`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- [[simulation.setup__ephemerides_time_seconds_flexible|_ephemerides_time_seconds_flexible]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1499-1499`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:49-49`
<!-- vulcan:connections:end -->

## Limitations
Resolution is truncated to 1 ms; microsecond epochs are rounded. Field values are converted with `Int(...)`, which throws `InexactError` for non-integral floats such as `year = 2024.5`. `DateTime` construction throws `ArgumentError` for out-of-range month or day. No time-zone handling exists; the result is interpreted as UTC by every caller.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 40.
