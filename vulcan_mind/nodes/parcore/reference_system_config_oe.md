---
id: parcore.reference_system_config_oe
label: OE
kind: struct
source:
  file: src/core/state/reference_system_config.jl
  symbol: OE
  lines:
  - 4
  - 12
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Base numeric types; the module is included directly by SimulationModel
    before the submodule block.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: elements
  type: OE
  units: m, -, rad
  description: Mutable named orbital-element record with semi-major axis, eccentricity,
    inclination, RAAN, argument of periapsis, true anomaly and mass.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# OE

## Purpose
`OE` is the named record for classical orbital elements inside the `ReferenceSystems` module. It gives configuration and reporting code a field-addressable form of the same quantities that `rvtoorbitalelement` returns as a static vector, so that scripts can write `oe.a` rather than index into position five of an anonymous tuple.

## Model & Assumptions
Fields are plain `Float64` with the package convention that angles are radians and lengths are metres. The record carries mass `m` alongside the six geometric elements, mirroring the seven-element static vector used in the aerobraking path. Being mutable, it is an in-place state holder rather than a value type, so two references to the same record observe each other's edits.

## Design & Implementation
The file is a small container module declaring six sibling records: `OE`, `cartesian`, `R_RA_DEC`, `H_LAN_LON`, `uDuNuE` and `clock`, all exported. Each is a `mutable struct` of `Float64` fields except `clock`, whose calendar fields are `Int64` with a fractional `second`. The module is included by `simulation_model.jl` in the utility block ahead of the submodule sequence, because configuration parsing and the reference-frame interface both need these names before the physics modules load.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Base numeric types; the module is included directly by SimulationModel before the submodule block. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `elements` | OE | m, -, rad | — | Mutable named orbital-element record with semi-major axis, eccentricity, inclination, RAAN, argument of periapsis, true anomaly and mass. |
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
No constructor validates ranges, so an eccentricity above one or a negative semi-major axis is accepted and only fails later in the conversion routines. Because the fields are untagged `Float64`, there is no compile-time protection against passing degrees where radians are expected. Mutability makes the record unsuitable for use as a dictionary key and defeats some of the escape analysis that the immutable static-vector form enjoys inside the integrator.

## Provenance
Mapped from `src/core/state/reference_system_config.jl:4-12`.
