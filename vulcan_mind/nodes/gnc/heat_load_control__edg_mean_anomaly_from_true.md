---
id: gnc.heat_load_control__edg_mean_anomaly_from_true
label: _edg_mean_anomaly_from_true
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_mean_anomaly_from_true
  lines:
  - 82
  - 82
inputs:
- id: nu
  type: Float64
  units: n/a
  required: true
  description: Positional argument `nu`.
- id: e
  type: Float64
  units: n/a
  required: true
  description: Positional argument `e`.
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
  description: Return value of `_edg_mean_anomaly_from_true`.
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

# _edg_mean_anomaly_from_true

## Purpose
Converts true anomaly to mean anomaly via Kepler's equation, giving a linear-in-time angle for computing elapsed time between two points on an elliptical orbit.

## Theory & Math
$$M = \left(E - e\sin E\right) \bmod 2\pi$$ with $E$ the eccentric anomaly (rad) and $e$ the eccentricity.

## Design & Implementation
`@inline` function `(nu::Float64, e::Float64)::Float64`. Calls `_edg_eccentric_anomaly_from_true(nu, e)` to obtain `E`, then returns `mod(E - e * sin(E), 2pi)`. Used twice by `_edg_drag_passage_duration` for the current and exit true anomalies.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nu` | Float64 | n/a | yes | Positional argument `nu`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_mean_anomaly_from_true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:96-96`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_eccentric_anomaly_from_true|_edg_eccentric_anomaly_from_true]] · `callers` · call · `src/gnc/control/heat_load_control.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
Inherits the elliptical-only validity of the eccentric-anomaly conversion. The `mod 2pi` wrap means differences between two mean anomalies must themselves be wrapped by the caller, which `_edg_drag_passage_duration` does.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 82.
