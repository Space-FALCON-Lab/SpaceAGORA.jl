---
id: core.reference_system_config_r_ra_dec
label: R_RA_DEC
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: R_RA_DEC
  lines:
  - 20
  - 20
inputs:
- id: r
  type: Float64
  units: n/a
  required: true
  description: Field `r`.
- id: RA
  type: Float64
  units: n/a
  required: true
  description: Field `RA`.
- id: dec
  type: Float64
  units: n/a
  required: true
  description: Field `dec`.
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
  type: R_RA_DEC
  units: n/a
  description: Constructed `R_RA_DEC`.
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

# R_RA_DEC

## Purpose

`R_RA_DEC` is the mutable spherical-coordinate container for a position expressed as range `r`, right ascension `RA` and declination `dec`, all `Float64`. It is the celestial-sphere counterpart to `cartesian` for pointing vectors and line-of-sight directions.

## Design & Implementation

It is a bare `mutable struct` with three fields, no inner constructor and no normalisation, so instances are made with the default positional constructor `R_RA_DEC(r, RA, dec)` and mutated field by field. Conversion to and from `cartesian` is not defined here; the module holds only the data shape, leaving the trigonometry to the astrodynamics code that consumes these types.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | Float64 | n/a | yes | Field `r`. |
| in | `RA` | Float64 | n/a | yes | Field `RA`. |
| in | `dec` | Float64 | n/a | yes | Field `dec`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | R_RA_DEC | n/a | — | Constructed `R_RA_DEC`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/reference_system_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Angle units are not encoded, so degrees versus radians and hours-of-right-ascension versus degrees are conventions the caller must track. Nothing constrains `dec` to the range $[-90°, 90°]$ or wraps `RA` into $[0°, 360°)$, and a range `r` of zero leaves the two angles undefined without any error being raised. The type also records no equinox or epoch, so J2000 and true-of-date directions share one representation.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 20.
