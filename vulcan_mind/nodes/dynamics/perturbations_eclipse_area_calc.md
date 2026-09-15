---
id: dynamics.perturbations_eclipse_area_calc
label: eclipse_area_calc
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: eclipse_area_calc
  lines:
  - 2277
  - 2277
inputs:
- id: r_sat
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_sat`.
- id: r_sun
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_sun`.
- id: rp
  type: Float64
  units: n/a
  required: true
  description: Positional argument `rp`.
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
  description: Return value of `eclipse_area_calc`. Returns `1.0` or `0.0` or `1.0
    - b^2 / a^2` or `1 - A / (π * a^2)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# eclipse_area_calc

## Purpose
Computes the fraction of the solar disc visible from the spacecraft, accounting for umbra, penumbra and annular geometry, following the Basilisk conical shadow model.

## Design & Implementation
Returns one if either vector is degenerate or the spacecraft is on the sunward side. It computes the apparent angular radii of Sun and planet and their apparent separation, returning zero for total eclipse, `1 - b²/a²` for annular, and for partial eclipse the fraction of the solar disc outside the planet's disc via the two-circle intersection area. Uses a Sun radius of 695,000 km.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_sat` | SVector{3, Float64} | n/a | yes | Positional argument `r_sat`. |
| in | `r_sun` | SVector{3, Float64} | n/a | yes | Positional argument `r_sun`. |
| in | `rp` | Float64 | n/a | yes | Positional argument `rp`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `eclipse_area_calc`. Returns `1.0` or `0.0` or `1.0 - b^2 / a^2` or `1 - A / (π * a^2)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_srp_cannonball_accel|srp_cannonball_accel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1137-1137`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Spherical planet and no atmospheric refraction; the partial-eclipse intersection formula assumes small angles and loses accuracy for very low orbits.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2277.
