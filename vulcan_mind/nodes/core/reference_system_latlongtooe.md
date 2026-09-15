---
id: core.reference_system_latlongtooe
label: latlongtoOE
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: latlongtoOE
  lines:
  - 286
  - 286
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
- id: gamma
  type: Any
  units: n/a
  required: true
  description: Positional argument `γ`.
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `α`.
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
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
  description: Return value of `latlongtoOE`. Returns `OE`.
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

# latlongtoOE

## Purpose
Converts a geodetic position plus a local velocity given by flight-path angle, azimuth and speed into classical orbital elements, for entry-style initial conditions.

## Design & Implementation
Builds the ellipsoidal Cartesian position, converts it to inertial with the legacy `r_pintor_i`, forms local zenith, east and north unit vectors in the planet-fixed frame, decomposes the velocity into north, east and up components with `γ` and `α`, rotates the resulting planet-fixed velocity to J2000 through `L_PI'`, and calls `rvtoorbitalelement`. Returns the first six elements.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `LATLONGH` | Any | n/a | yes | Positional argument `LATLONGH`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `gamma` | Any | n/a | yes | Positional argument `γ`. |
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `latlongtoOE`. Returns `OE`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/core/interfaces/reference_system.jl:306-306`
- `callees` → [[core.reference_system_r_pintor_i|r_pintor_i]] · `callers` · call · `src/core/interfaces/reference_system.jl:324-324`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/core/interfaces/reference_system.jl:347-347`
<!-- vulcan:connections:end -->

## Limitations
It prints every intermediate value with `println`, which is debugging output left in production code; the planet-fixed velocity is rotated without the transport term, so the resulting inertial velocity omits the planet's rotation and the elements are for a non-rotating planet.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 286.
