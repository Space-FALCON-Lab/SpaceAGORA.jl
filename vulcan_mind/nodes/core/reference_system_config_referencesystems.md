---
id: core.reference_system_config_referencesystems
label: ReferenceSystems
kind: module
source:
  file: src/core/state/reference_system_config.jl
  symbol: ReferenceSystems
  lines:
  - 1
  - 1
inputs:
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
  description: Value produced by this symbol.
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

# ReferenceSystems

## Purpose

`ReferenceSystems` is a small module that declares the coordinate and time container types used across SpaceAGORA for expressing spacecraft state in different reference systems. It exports `OE`, `cartesian`, `R_RA_DEC`, `H_LAN_LON`, `uDuNuE` and `clock`.

## Design & Implementation

The module body is nothing but six `mutable struct` declarations and a single `export` line; there are no constructors, conversion methods or validation. Every field is a concrete `Float64` except `clock`, whose calendar fields are `Int64` with a `Float64` seconds field. Making the structs mutable lets propagation and conversion code update a state container in place instead of allocating a replacement each step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
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

Because the module defines only data and no conversion functions, the transformations between these representations live elsewhere and there is nothing here to keep them consistent. The generic names `cartesian` and `clock` are lower case and unqualified, so `using ReferenceSystems` risks colliding with identically named bindings from other packages. No units, epochs or frames are recorded on the types, so an ECI position and an ECEF position are indistinguishable at the type level.

## Provenance
Mapped from `src/core/state/reference_system_config.jl` line 1.
