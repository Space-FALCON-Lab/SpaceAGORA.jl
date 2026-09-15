---
id: gnc.targeting_control__edg_orbit_metrics_from_rv
label: _edg_orbit_metrics_from_rv
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_orbit_metrics_from_rv
  lines:
  - 292
  - 292
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: Tuple
  units: n/a
  description: Return value of `_edg_orbit_metrics_from_rv`. Returns `(energy=energy,
    periapsis=NaN, apoapsis=Inf)` or `(energy=energy, periapsis=NaN, apoapsis=NaN)`
    or `(energy=energy, periapsis=a * (1.0 - e), apoapsis=a * (1.0 + e))`.
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

# _edg_orbit_metrics_from_rv

## Purpose
Summarises an inertial state as specific orbital energy and periapsis and apoapsis radii, the quantities targeting compares against.

## Theory & Math
$$
\epsilon = \tfrac{1}{2}v^2 - \frac{\mu}{r},\qquad r_p = a(1 - e),\qquad r_a = a(1 + e)
$$

## Design & Implementation
Computes energy as `v²/2 - μ/r`. If it is non-negative or non-finite the orbit is unbound and the tuple carries `NaN` periapsis and infinite apoapsis. Otherwise it converts to elements and returns `a(1-e)` and `a(1+e)`, substituting `NaN` for both if the elements are non-finite.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_edg_orbit_metrics_from_rv`. Returns `(energy=energy, periapsis=NaN, apoapsis=Inf)` or `(energy=energy, periapsis=NaN, apoapsis=NaN)` or `(energy=energy, periapsis=a * (1.0 - e), apoapsis=a * (1.0 + e))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:756-756`
- [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:699-699`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/targeting_control.jl:299-299`
<!-- vulcan:connections:end -->

## Limitations
For an orbit that is barely bound the element conversion is ill-conditioned, so apoapsis can swing by large amounts between neighbouring switch candidates, which the certification step has to absorb.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 292.
