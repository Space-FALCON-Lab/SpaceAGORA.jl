---
id: core.runtime_types_initialparameters
label: InitialParameters
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: InitialParameters
  lines:
  - 41
  - 41
inputs:
- id: M
  type: Mission
  units: n/a
  required: false
  description: Field `M` (default `Mission()`).
- id: gm
  type: LegacyGravityModelCode
  units: n/a
  required: false
  description: Field `gm` (default `LegacyGravityConstant`).
- id: dm
  type: LegacyDensityModelCode
  units: n/a
  required: false
  description: Field `dm` (default `LegacyDensityConstant`).
- id: wm
  type: Int64
  units: n/a
  required: false
  description: Field `wm` (default `0`).
- id: am
  type: LegacyAerodynamicModelCode
  units: n/a
  required: false
  description: Field `am` (default `LegacyAeroCdClConstant`).
- id: tm
  type: LegacyThermalModelCode
  units: n/a
  required: false
  description: Field `tm` (default `LegacyThermalConvectiveRadiative`).
- id: cm
  type: Int64
  units: n/a
  required: false
  description: Field `cm` (default `0`).
- id: tc
  type: LegacyThrustControlCode
  units: n/a
  required: false
  description: Field `tc` (default `LegacyThrustNone`).
- id: mc
  type: Int64
  units: n/a
  required: false
  description: Field `mc` (default `0`).
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
  type: InitialParameters
  units: n/a
  description: Constructed `InitialParameters` (keyword constructor via @kwdef).
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

# InitialParameters

## Purpose
Legacy model-selection record mapping integer or enum codes to gravity, density, aerodynamic, thermal, and thrust-control model choices, with a compatibility constructor for integer-coded callers.

## Design & Implementation
`@kwdef struct InitialParameters` with `M::Mission`, enum fields `gm::LegacyGravityModelCode` (default `LegacyGravityConstant`), `dm::LegacyDensityModelCode` (`LegacyDensityConstant`), `am::LegacyAerodynamicModelCode` (`LegacyAeroCdClConstant`), `tm::LegacyThermalModelCode` (`LegacyThermalConvectiveRadiative`), `tc::LegacyThrustControlCode` (`LegacyThrustNone`), and integer fields `wm`, `cm`, `mc`. A nine-positional-argument outer constructor accepts `Union{Enum, Integer}` for each coded slot and routes them through `_compat_enum_parse` before calling the keyword constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `M` | Mission | n/a | no | Field `M` (default `Mission()`). |
| in | `gm` | LegacyGravityModelCode | n/a | no | Field `gm` (default `LegacyGravityConstant`). |
| in | `dm` | LegacyDensityModelCode | n/a | no | Field `dm` (default `LegacyDensityConstant`). |
| in | `wm` | Int64 | n/a | no | Field `wm` (default `0`). |
| in | `am` | LegacyAerodynamicModelCode | n/a | no | Field `am` (default `LegacyAeroCdClConstant`). |
| in | `tm` | LegacyThermalModelCode | n/a | no | Field `tm` (default `LegacyThermalConvectiveRadiative`). |
| in | `cm` | Int64 | n/a | no | Field `cm` (default `0`). |
| in | `tc` | LegacyThrustControlCode | n/a | no | Field `tc` (default `LegacyThrustNone`). |
| in | `mc` | Int64 | n/a | no | Field `mc` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | InitialParameters | n/a | — | Constructed `InitialParameters` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- `callees` → [[core.compat_model_codes__compat_enum_parse|_compat_enum_parse]] · `callers` · call · `src/core/types/runtime_types.jl:67-67`
- `callees` → [[core.runtime_types_mission|Mission]] · `callers` · call · `src/core/types/runtime_types.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
The positional constructor throws whatever `_compat_enum_parse` raises for an unknown integer code, with no context about which field failed. `wm`, `cm`, and `mc` remain untyped integers with undocumented meanings. This record is only consumed by the legacy Python-port paths.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 41.
