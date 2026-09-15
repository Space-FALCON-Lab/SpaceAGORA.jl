---
id: gnc.targeting_control__edg_ephemeris_time
label: _edg_ephemeris_time
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_ephemeris_time
  lines:
  - 315
  - 315
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
  type: Float64
  units: n/a
  description: Return value of `_edg_ephemeris_time`.
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

# _edg_ephemeris_time

## Purpose
Converts simulation elapsed time to SPICE ephemeris time when the run carries an epoch, falling back to elapsed time in harnesses that do not.

## Design & Implementation
Returns `p.shared_buffers.et_start[] + t_abs` when both properties exist, otherwise `t_abs`. `@inline`. The property checks let the controller run in a minimal test `ODEParams` without shared buffers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t_abs` | Float64 | n/a | yes | Positional argument `t_abs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_ephemeris_time`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:69-69`
- [[gnc.targeting_control__edg_planet_frame_lpi|_edg_planet_frame_lpi]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:325-325`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:331-331`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The fallback silently treats elapsed seconds as ephemeris time, so in a harness without `et_start` any SPICE-based frame is evaluated at the wrong epoch.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 315.
