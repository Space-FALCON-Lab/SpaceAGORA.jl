---
id: simulation.setup__ephemerides_time_seconds_flexible
label: _ephemerides_time_seconds_flexible
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemerides_time_seconds_flexible
  lines:
  - 1491
  - 1491
inputs:
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: ephemerides_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `ephemerides_model`.
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
  description: Return value of `_ephemerides_time_seconds_flexible`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _ephemerides_time_seconds_flexible

## Purpose
Computes the ephemeris start time for a configuration even when the ephemerides model lacks the standard `ephemerides_time_seconds` method, as older or partially loaded models do.

## Design & Implementation
Uses `applicable` to prefer the generic method. Otherwise it branches on the model type name: `SpiceEphemeridesModel` converts the initial time to UTC and calls `utc2et` under the SPICE lock; `SimpleEphemeridesModel` adds the elapsed milliseconds since the J2000 UTC epoch to the model's reference seconds. Any other type raises `MethodError`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `ephemerides_model` | Any | n/a | yes | Positional argument `ephemerides_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_ephemerides_time_seconds_flexible`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1712-1712`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:1507-1507`
- `callees` → [[environment.simple_ephemerides__initial_time_datetime|_initial_time_datetime]] · `callers` · call · `src/simulation/engine/setup.jl:1499-1499`
- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/engine/setup.jl:1493-1493`
<!-- vulcan:connections:end -->

## Limitations
Dispatching on `nameof(typeof(...))` rather than the type itself is a workaround for load-order issues and would misfire on a user type that happened to share a name.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1491.
