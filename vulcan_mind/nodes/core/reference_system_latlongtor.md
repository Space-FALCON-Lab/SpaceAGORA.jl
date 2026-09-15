---
id: core.reference_system_latlongtor
label: latlongtor
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: latlongtor
  lines:
  - 267
  - 267
inputs:
- id: LATLONGH
  type: Any
  units: n/a
  required: true
  description: Positional argument `LATLONGH`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: alpha_g0
  type: Any
  units: n/a
  required: true
  description: Positional argument `α_g0`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: t0
  type: Any
  units: n/a
  required: true
  description: Positional argument `t0`.
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
  type: AbstractArray
  units: n/a
  description: Return value of `latlongtor`. Returns `[x, y, z]`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# latlongtor

## Purpose
Converts geodetic latitude, longitude and altitude to an inertial position at a given time, accounting for the planet's rotation since a reference epoch.

## Theory & Math
$$
N = \frac{a}{\sqrt{1 - e^2 \sin^2\phi}},\qquad x = (N + h)\cos\phi\cos\alpha,\quad y = (N + h)\cos\phi\sin\alpha,\quad z = \left((1 - e^2) N + h\right)\sin\phi
$$

## Design & Implementation
Computes the first eccentricity from the equatorial and polar radii, rotates the longitude by `α_g0 + ω_z (t - t0)` to obtain inertial right ascension, evaluates the prime-vertical radius `N`, and returns the ellipsoidal Cartesian coordinates as a `Vector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `LATLONGH` | Any | n/a | yes | Positional argument `LATLONGH`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `alpha_g0` | Any | n/a | yes | Positional argument `α_g0`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `t0` | Any | n/a | yes | Positional argument `t0`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `latlongtor`. Returns `[x, y, z]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the z component of the spin vector is used, so it is correct only for a pole-aligned spin model; the returned vector is heap-allocated.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 267.
