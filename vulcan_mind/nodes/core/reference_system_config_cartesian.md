---
id: core.reference_system_config_cartesian
label: cartesian
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: cartesian
  lines:
  - 14
  - 14
inputs:
- id: x
  type: Float64
  units: n/a
  required: true
  description: Field `x`.
- id: y
  type: Float64
  units: n/a
  required: true
  description: Field `y`.
- id: z
  type: Float64
  units: n/a
  required: true
  description: Field `z`.
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
  type: cartesian
  units: n/a
  description: Constructed `cartesian`.
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

# cartesian

## Purpose

`cartesian` is the mutable three-component container for a rectangular coordinate triple, holding `x`, `y` and `z` as `Float64`. It is the plain position or vector representation that the other reference-system types in this module are converted to and from.

## Design & Implementation

The declaration is a bare `mutable struct` with three `Float64` fields and no inner constructor, so the default positional constructor `cartesian(x, y, z)` is the only way to build one. Mutability means propagation code can assign `c.x = ...` in place rather than allocating a fresh instance per step. The type is exported from `ReferenceSystems` alongside `OE`, `R_RA_DEC`, `H_LAN_LON`, `uDuNuE` and `clock`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Float64 | n/a | yes | Field `x`. |
| in | `y` | Float64 | n/a | yes | Field `y`. |
| in | `z` | Float64 | n/a | yes | Field `z`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | cartesian | n/a | — | Constructed `cartesian`. |
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

The struct carries no frame tag and no unit annotation, so an ECI position in kilometres and a body-frame offset in metres share the same type and can be mixed without any error. There is no arithmetic defined on it, so callers must unpack the fields before doing vector maths. The lower-case name is unusual for a Julia type and is easy to confuse with a function when it appears in exported scope.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 14.
