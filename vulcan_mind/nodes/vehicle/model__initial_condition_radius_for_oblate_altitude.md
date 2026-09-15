---
id: vehicle.model__initial_condition_radius_for_oblate_altitude
label: _initial_condition_radius_for_oblate_altitude
kind: function
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: _initial_condition_radius_for_oblate_altitude
  lines:
  - 146
  - 146
inputs:
- id: target_altitude
  type: Float64
  units: n/a
  required: true
  description: Positional argument `target_altitude`.
- id: u_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `u_pp`.
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
  type: Float64
  units: n/a
  description: Return value of `_initial_condition_radius_for_oblate_altitude`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# _initial_condition_radius_for_oblate_altitude

## Purpose
Finds the geocentric radius along a planet-fixed direction at which the geodetic altitude equals `target_altitude`, so apoapsis and periapsis altitudes specified above the ellipsoid can be turned into orbital radii for `InitialCondition(ra=, rp=)`.

## Design & Implementation
Takes `target_altitude::Float64` (m), `u_pp::SVector{3,Float64}` and `planet`. Negative targets throw `ArgumentError("Oblate InitialCondition altitudes must be nonnegative; got ... m.")`. The bracket starts at `lo = _initial_condition_oblate_surface_radius(u_pp, planet)` and `hi = lo + target + |Rp_e-Rp_p| + 1`, expanding `hi` by `max(target, |Rp_e-Rp_p|, 1)` until the altitude at `hi` reaches the target. Exactly 80 bisection iterations then narrow `[lo, hi]` using `_initial_condition_oblate_altitude`, and the midpoint is returned in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `target_altitude` | Float64 | n/a | yes | Positional argument `target_altitude`. |
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_initial_condition_radius_for_oblate_altitude`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/vehicle/spacecraft/model.jl:198-198`
- `callees` → [[vehicle.model__initial_condition_apsis_direction_ii|_initial_condition_apsis_direction_ii]] · `callers` · call · `src/vehicle/spacecraft/model.jl:203-203`
- `callees` → [[vehicle.model__initial_condition_lpi|_initial_condition_lpi]] · `callers` · call · `src/vehicle/spacecraft/model.jl:201-201`
- `callees` → [[vehicle.model__initial_condition_oblate_altitude|_initial_condition_oblate_altitude]] · `callers` · call · `src/vehicle/spacecraft/model.jl:156-156`
- `callees` → [[vehicle.model__initial_condition_oblate_surface_radius|_initial_condition_oblate_surface_radius]] · `callers` · call · `src/vehicle/spacecraft/model.jl:154-154`
- `callees` → [[vehicle.model_initialcondition|InitialCondition]] · `callers` · call · `src/vehicle/spacecraft/model.jl:172-172`
<!-- vulcan:connections:end -->

## Limitations
Inherits the altitude formula discrepancy of `_initial_condition_oblate_altitude`, so the computed radius corresponds to that formula's altitude rather than the textbook geodetic altitude. The bisection count is fixed at 80 regardless of bracket width, and the expansion loop has no upper bound; a NaN altitude (degenerate direction) makes the comparison false and returns the initial bracket midpoint silently.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 146.
