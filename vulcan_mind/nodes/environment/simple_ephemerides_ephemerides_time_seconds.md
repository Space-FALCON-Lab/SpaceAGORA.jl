---
id: environment.simple_ephemerides_ephemerides_time_seconds
label: ephemerides_time_seconds
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: ephemerides_time_seconds
  lines:
  - 52
  - 52
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
  type: Float64
  units: n/a
  description: Return value of `ephemerides_time_seconds`.
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

# ephemerides_time_seconds

## Purpose
Converts the mission's `initial_time` into the scalar ephemeris time (seconds) that the selected ephemeris model uses as its timeline origin: SPICE ephemeris time (TDB seconds past J2000) for `SpiceEphemeridesModel`, or UTC seconds past J2000 offset by `reference_epoch_seconds` for `SimpleEphemeridesModel`.

## Design & Implementation
Both methods first call `_initial_time_datetime(initial_time)`. The SPICE method wraps it with `from_utc`, then under `lock(SPICE_LOCK)` calls `utc2et(to_utc(start_epoch))`, which applies leap seconds and the TDB offset from the loaded LSK kernel. The simple method computes `Dates.value(start_time - _J2000_UTC)` in milliseconds (with `_J2000_UTC = DateTime(2000,1,1,12,0,0)`) and returns `model.reference_epoch_seconds + elapsed_ms / 1000.0`. Return type is `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `ephemerides_time_seconds`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:510-510`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/setup.jl:1493-1493`
- [[simulation.planet_frame_init_affect_bang|init_affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- [[simulation.setup__ephemerides_time_seconds_flexible|_ephemerides_time_seconds_flexible]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1493-1493`
- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:219-219`
- [[vehicle.model__initial_condition_lpi|_initial_condition_lpi]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:99-99`

**Downstream**

- `callees` → [[environment.simple_ephemerides__initial_time_datetime|_initial_time_datetime]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:53-53`
<!-- vulcan:connections:end -->

## Limitations
The simple timeline ignores leap seconds and TDB-UTC, so it differs from SPICE ET by roughly 64-69 s for modern epochs; the two models therefore are not interchangeable at the second level. The SPICE method throws if no leap-seconds kernel is furnished. Millisecond resolution of `DateTime` limits both paths. Locking on every call is negligible here since it is invoked once per simulation setup.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 52.
