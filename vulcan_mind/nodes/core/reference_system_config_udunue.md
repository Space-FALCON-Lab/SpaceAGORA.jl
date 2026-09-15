---
id: core.reference_system_config_udunue
label: uDuNuE
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: uDuNuE
  lines:
  - 32
  - 32
inputs:
- id: uD
  type: Float64
  units: n/a
  required: true
  description: Field `uD`.
- id: uN
  type: Float64
  units: n/a
  required: true
  description: Field `uN`.
- id: uE
  type: Float64
  units: n/a
  required: true
  description: Field `uE`.
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
  type: uDuNuE
  units: n/a
  description: Constructed `uDuNuE`.
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

# uDuNuE

## Purpose

`uDuNuE` is the mutable container for a vector resolved onto the local topocentric Down-North-East unit triad, storing the three components `uD`, `uN` and `uE` as `Float64`. It expresses a direction or velocity relative to the local horizon at a point on the reference body.

## Design & Implementation

Like the other types in `ReferenceSystems` it is a bare `mutable struct` with no inner constructor, no validation and no arithmetic, constructed positionally as `uDuNuE(uD, uN, uE)` and updated in place. The Down-North-East ordering is a NED-family convention with Down first, so callers converting from a standard NED vector must reorder the components rather than copying them straight across.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `uD` | Float64 | n/a | yes | Field `uD`. |
| in | `uN` | Float64 | n/a | yes | Field `uN`. |
| in | `uE` | Float64 | n/a | yes | Field `uE`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | uDuNuE | n/a | — | Constructed `uDuNuE`. |
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

The triad it refers to is implicit: the type does not record the latitude and longitude of the origin point the frame is anchored at, so two instances anchored at different ground points are indistinguishable. No unit is attached, so it serves equally as a dimensionless unit vector or as a velocity in metres per second. Nothing normalises the triple, so a value nominally intended as a unit direction may not have unit norm.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 32.
