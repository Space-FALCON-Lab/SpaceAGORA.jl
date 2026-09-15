---
id: gnc.heat_load_control__edg_eccentric_anomaly_from_true
label: _edg_eccentric_anomaly_from_true
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_eccentric_anomaly_from_true
  lines:
  - 77
  - 77
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
  description: Return value of `_edg_eccentric_anomaly_from_true`.
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

# _edg_eccentric_anomaly_from_true

## Purpose
Converts true anomaly to eccentric anomaly for an elliptical orbit so the drag-passage duration can be found through Kepler's equation.

## Theory & Math
$$E = \operatorname{atan2}\!\left(\sqrt{1-e^2}\,\sin\nu,\; e + \cos\nu\right) \bmod 2\pi$$ where $\nu$ is the true anomaly (rad), $e$ the eccentricity, and $E$ the eccentric anomaly (rad).

## Design & Implementation
`@inline` function `(nu::Float64, e::Float64)::Float64`. Uses the two-argument arctangent `E = atan(sqrt(max(0, 1 - e^2)) * sin(nu), e + cos(nu))`, which is quadrant-safe, then returns `mod(E, 2pi)` in radians.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nu` | Float64 | n/a | yes | Positional argument `nu`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_eccentric_anomaly_from_true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_mean_anomaly_from_true|_edg_mean_anomaly_from_true]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:83-83`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Valid only for `0 <= e < 1`; the `max(0, 1 - e^2)` guard prevents a domain error for hyperbolic inputs but returns a meaningless angle. The caller must check eccentricity before use.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 77.
