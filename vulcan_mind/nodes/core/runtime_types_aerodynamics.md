---
id: core.runtime_types_aerodynamics
label: Aerodynamics
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Aerodynamics
  lines:
  - 119
  - 119
inputs:
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `0.0`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `0.0`).
- id: thermal_accomodation_factor
  type: Float64
  units: n/a
  required: false
  description: Field `thermal_accomodation_factor` (default `0.0`).
- id: reflection_coefficient
  type: Float64
  units: n/a
  required: false
  description: Field `reflection_coefficient` (default `0.0`).
- id: thermal_contact
  type: Float64
  units: n/a
  required: false
  description: Field `thermal_contact` (default `0.0`).
- id: heat_rate_limit
  type: Float64
  units: n/a
  required: false
  description: Field `heat_rate_limit` (default `0.0`).
- id: heat_load_limit
  type: Float64
  units: n/a
  required: false
  description: Field `heat_load_limit` (default `0.0`).
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
  type: Aerodynamics
  units: n/a
  description: Constructed `Aerodynamics` (keyword constructor via @kwdef).
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

# Aerodynamics

## Purpose
Legacy mutable record of scalar aerodynamic and thermal-limit parameters carried by the `Model` struct for single-vehicle aerobraking studies.

## Design & Implementation
`@kwdef mutable struct Aerodynamics` with seven `Float64` fields all defaulting to `0.0`: `δ` and `α` (angle-of-attack-related angles, rad), `thermal_accomodation_factor` (dimensionless, 0-1), `reflection_coefficient`, `thermal_contact`, `heat_rate_limit` (W/cm^2), and `heat_load_limit` (J/cm^2). Being mutable, it can be updated in place by legacy control code without rebuilding the enclosing `Model`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `delta` | Float64 | n/a | no | Field `δ` (default `0.0`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `0.0`). |
| in | `thermal_accomodation_factor` | Float64 | n/a | no | Field `thermal_accomodation_factor` (default `0.0`). |
| in | `reflection_coefficient` | Float64 | n/a | no | Field `reflection_coefficient` (default `0.0`). |
| in | `thermal_contact` | Float64 | n/a | no | Field `thermal_contact` (default `0.0`). |
| in | `heat_rate_limit` | Float64 | n/a | no | Field `heat_rate_limit` (default `0.0`). |
| in | `heat_load_limit` | Float64 | n/a | no | Field `heat_load_limit` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Aerodynamics | n/a | — | Constructed `Aerodynamics` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_model|Model]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:160-160`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Zero defaults are not physically meaningful (a zero accommodation factor implies no energy transfer), so every field must be set explicitly by the caller. The struct is untyped in its units and has no validation; negative limits are accepted. Modern multibody paths use per-link coefficients from `SpacecraftModel` instead and largely ignore this record.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 119.
