---
id: gnc.heat_load_control__edg_heat_load_track_env
label: _edg_heat_load_track_env
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_track_env
  lines:
  - 316
  - 316
inputs:
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: j
  type: Int
  units: n/a
  required: true
  description: Positional argument `j`.
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
  description: Return value of `_edg_heat_load_track_env`. Returns `(`.
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

# _edg_heat_load_track_env

## Purpose
Packs the atmospheric and kinematic state at a single track node into the NamedTuple layout consumed by the heat-rate and structural-load root solvers.

## Design & Implementation
Takes `track` and an index `j::Int`. It clamps `speed = max(track.speed[j], 0)` and `rho = max(track.rho[j], 0)`, then returns `(altitude_m = track.h[j], rho, temperature = track.temperature[j], speed, molecular_speed_ratio = track.speed_ratio[j], dynamic_pressure = 0.5 rho speed^2)` with dynamic pressure in Pa.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `j` | Int | n/a | yes | Positional argument `j`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_heat_load_track_env`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:353-353`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No bounds check on `j`; an index outside the track raises `BoundsError`. Temperature and speed ratio are passed through unclamped even though negative values would be unphysical. Allocates a small tuple per call inside the per-node constraint loop.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 316.
