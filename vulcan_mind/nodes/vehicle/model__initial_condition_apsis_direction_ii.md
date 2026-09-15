---
id: vehicle.model__initial_condition_apsis_direction_ii
label: _initial_condition_apsis_direction_ii
kind: function
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: _initial_condition_apsis_direction_ii
  lines:
  - 110
  - 110
inputs:
- id: i
  type: Float64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: omega
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ω`.
- id: Omega
  type: Float64
  units: n/a
  required: true
  description: Positional argument `Ω`.
- id: nu
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ν`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_initial_condition_apsis_direction_ii`.
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

# _initial_condition_apsis_direction_ii

## Purpose
Computes the inertial unit vector pointing from the planet centre toward a point on the orbit at true anomaly `ν`, used by the oblate `InitialCondition(planet; ra, hp, ...)` constructor to find the apoapsis (`ν = π`) and periapsis (`ν = 0`) directions before converting altitudes to radii.

## Theory & Math
$\hat{u}_{II} = Q^{T}\begin{pmatrix}\cos\nu\\ \sin\nu\\ 0\end{pmatrix}/\|\cdot\|$ where $Q = R_3(\omega)R_1(i)R_3(\Omega)$ is the inertial-to-perifocal rotation, $\Omega$ the RAAN, $i$ inclination, $\omega$ argument of periapsis and $\nu$ true anomaly (all rad).

## Design & Implementation
`@inline` function of `i`, `ω`, `Ω`, `ν` (all `Float64`, radians) returning `SVector{3,Float64}`. It forms the perifocal direction `q = (cos ν, sin ν, 0)` and the 3-by-3 rotation `Q` from inertial to perifocal coordinates (rows built from the classical 3-1-3 sequence in `Ω`, `i`, `ω`), then returns `Q' * q` normalised by its norm. Because `Q` is orthonormal the normalisation is a numerical safeguard only.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `i` | Float64 | n/a | yes | Positional argument `i`. |
| in | `omega` | Float64 | n/a | yes | Positional argument `ω`. |
| in | `Omega` | Float64 | n/a | yes | Positional argument `Ω`. |
| in | `nu` | Float64 | n/a | yes | Positional argument `ν`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_initial_condition_apsis_direction_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:203-203`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The matrix `Q` is written out element by element, so a transcription error would not be caught by construction; it matches the standard perifocal-to-inertial transpose. Angles must already be in radians; the caller converts from degrees. The function is only exercised for `ν` in {0, π} by the constructor, so other anomalies are untested here.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 110.
