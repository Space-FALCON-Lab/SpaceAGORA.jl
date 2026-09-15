---
id: gnc.targeting_control__edg_planet_frame_lpi
label: _edg_planet_frame_lpi
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_planet_frame_lpi
  lines:
  - 322
  - 322
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t_abs
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t_abs`.
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
  description: Return value of `_edg_planet_frame_lpi`. Returns `planet_frame_lpi(planet,
    _edg_ephemeris_time(p, t_abs), ephemerides_model)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _edg_planet_frame_lpi

## Purpose
Fetches the inertial-to-planet-fixed rotation for the controller's prediction at a given elapsed time.

## Design & Implementation
Reads the planet and ephemerides model from the configuration and calls `planet_frame_lpi` at the converted ephemeris time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t_abs` | Float64 | n/a | yes | Positional argument `t_abs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_planet_frame_lpi`. Returns `planet_frame_lpi(planet, _edg_ephemeris_time(p, t_abs), ephemerides_model)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:332-332`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/gnc/control/targeting_control.jl:325-325`
- `callees` → [[gnc.targeting_control__edg_ephemeris_time|_edg_ephemeris_time]] · `callers` · call · `src/gnc/control/targeting_control.jl:325-325`
<!-- vulcan:connections:end -->

## Limitations
Called once per prediction step of the internal RK4 integration, so a SPICE-backed ephemerides model makes the switch solve dominated by kernel lookups.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 322.
