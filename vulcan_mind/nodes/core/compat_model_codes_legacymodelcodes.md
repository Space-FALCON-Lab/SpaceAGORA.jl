---
id: core.compat_model_codes_legacymodelcodes
label: LegacyModelCodes
kind: module
source:
  file: src/core/types/compat_model_codes.jl
  symbol: LegacyModelCodes
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

# LegacyModelCodes

## Purpose

`LegacyModelCodes` defines the `Int8`-backed enumerations that name the physics model selections inherited from the pre-Julia aerobraking tooling: gravity, atmospheric density, aerodynamics, thermal and thrust control. It exists so those fields can be strongly typed while legacy input decks that store bare integers keep working.

## Design & Implementation

Five `@enum ... ::Int8` declarations fix the wire values: gravity 0-3 (constant, inverse-square, inverse-square plus J2, GRAM), density 0-4 (constant, exponential, none, GRAM, NRLMSISE), aerodynamics 0-2 (constant Cd/Cl, diffusive, no ballistic axial), thermal 1-2 (convective-radiative, Maxwellian) and thrust control 0-2 (none, aerobraking maneuver, drag-passage firing). The module exports every code type and member, provides `_compat_enum_parse` for integer-to-enum conversion, and overloads `Base.:(==)` in both argument orders between each enum and `Integer` so expressions such as `ip.gm == 1` still evaluate correctly.

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

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/compat_model_codes.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The numeric values are frozen by the legacy file format and cannot be reordered without breaking existing input decks. Thermal codes deliberately start at 1, so 0 is invalid there while it is valid everywhere else. The equality overloads compare only against `Integer`, and because they cross enum families a gravity code and a density code with the same underlying value both compare equal to the same integer, so the type must be checked separately.

## Provenance
Mapped from `src/core/types/compat_model_codes.jl` line 1.
