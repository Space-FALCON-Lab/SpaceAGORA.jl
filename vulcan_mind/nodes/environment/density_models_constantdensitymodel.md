---
id: environment.density_models_constantdensitymodel
label: ConstantDensityModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: ConstantDensityModel
  lines:
  - 170
  - 170
inputs:
- id: density_kg_m3
  type: Float64
  units: n/a
  required: true
  description: Field `density_kg_m3`.
- id: temperature_k
  type: Float64
  units: n/a
  required: true
  description: Field `temperature_k`.
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
  type: ConstantDensityModel
  units: n/a
  description: Constructed `ConstantDensityModel` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# ConstantDensityModel

## Purpose
A uniform density and temperature everywhere, for unit tests and analytic checks where the atmosphere must be trivially known.

## Design & Implementation
A `@kwdef struct` with `density_kg_m3` and `temperature_k`. Its `getDensity` ignores position and time and returns the two constants with zero wind.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_kg_m3` | Float64 | n/a | yes | Field `density_kg_m3`. |
| in | `temperature_k` | Float64 | n/a | yes | Field `temperature_k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstantDensityModel | n/a | — | Constructed `ConstantDensityModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It has no `getDensityBatch!` specialisation, so batch callers fall through to the generic per-element path; and it does not declare itself as vanishing above the entry interface, so the split solver keeps aerodynamics engaged on coast arcs.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 170.
