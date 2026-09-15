---
id: core.reference_system_config_h_lan_lon
label: H_LAN_LON
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: H_LAN_LON
  lines:
  - 26
  - 26
inputs:
- id: h
  type: Float64
  units: n/a
  required: true
  description: Field `h`.
- id: LAT
  type: Float64
  units: n/a
  required: true
  description: Field `LAT`.
- id: LON
  type: Float64
  units: n/a
  required: true
  description: Field `LON`.
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
  type: H_LAN_LON
  units: n/a
  description: Constructed `H_LAN_LON`.
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

# H_LAN_LON

## Purpose

`H_LAN_LON` is the mutable geodetic-style container holding altitude `h`, latitude `LAT` and longitude `LON`, all as `Float64`. It is the ground-track representation of a spacecraft position, used when a state has to be reported relative to the rotating body rather than in an inertial frame.

## Design & Implementation

The type is a bare `mutable struct` with three fields and no inner constructor or validation, built via the default positional constructor `H_LAN_LON(h, LAT, LON)`. Mutability allows in-place update as a trajectory is walked. Note that the field names disagree with the type name: the type says `LAN` while the second field is `LAT`, so the type name reads as an artefact rather than a description of its contents.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `h` | Float64 | n/a | yes | Field `h`. |
| in | `LAT` | Float64 | n/a | yes | Field `LAT`. |
| in | `LON` | Float64 | n/a | yes | Field `LON`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | H_LAN_LON | n/a | — | Constructed `H_LAN_LON`. |
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

No units are attached, so whether `LAT` and `LON` are degrees or radians, and whether `h` is metres or kilometres, is a convention held entirely by the calling code. Nothing wraps longitude into a canonical range or clamps latitude to the poles, so out-of-range values propagate silently. The type does not record which ellipsoid or reference body the altitude is measured against, so geodetic and geocentric altitudes are indistinguishable.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 26.
