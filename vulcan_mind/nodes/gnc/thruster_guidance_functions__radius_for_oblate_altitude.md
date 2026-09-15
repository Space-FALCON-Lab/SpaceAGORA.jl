---
id: gnc.thruster_guidance_functions__radius_for_oblate_altitude
label: _radius_for_oblate_altitude
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _radius_for_oblate_altitude
  lines:
  - 66
  - 66
inputs:
- id: target_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `target_altitude_m`.
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
  description: Return value of `_radius_for_oblate_altitude`.
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

# _radius_for_oblate_altitude

## Purpose
Inverts `_oblate_altitude_from_radius` numerically: given a target geodetic altitude along a planet-fixed direction, it finds the geocentric radius whose geodetic altitude equals that target, so the periapsis-raise burn can target an altitude above the oblate surface.

## Design & Implementation
Arguments are `target_altitude_m::Float64`, `u_pp::SVector{3,Float64}` and `planet`. Negative targets return `NaN`. The bracket starts at `lo = _oblate_surface_radius(u_pp, planet)` and `hi = lo + target + |Rp_e-Rp_p| + 1.0`, then `hi` is widened by `max(target, |Rp_e-Rp_p|, 1.0)` until `_oblate_altitude_from_radius(hi)` reaches the target. A fixed 80-iteration bisection follows, moving `lo` up when the midpoint altitude is below target and `hi` down otherwise, and the midpoint `0.5*(lo+hi)` is returned. Because altitude is monotone in radius along a fixed ray, the bracket is guaranteed to contain the root.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `target_altitude_m` | Float64 | n/a | yes | Positional argument `target_altitude_m`. |
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_radius_for_oblate_altitude`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:183-183`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- `callees` → [[gnc.thruster_guidance_functions__oblate_altitude_from_radius|_oblate_altitude_from_radius]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:70-70`
- `callees` → [[gnc.thruster_guidance_functions__oblate_surface_radius|_oblate_surface_radius]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
80 bisections on a bracket of a few kilometres gives far better than micrometre resolution, but the count is fixed rather than tolerance-driven, so it always costs 80 altitude evaluations. The widening loop has no iteration cap; a NaN altitude (from a degenerate `u_pp`) makes the comparison false and skips widening, then bisection converges to the initial bracket midpoint silently. `u_pp` is assumed normalised.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 66.
