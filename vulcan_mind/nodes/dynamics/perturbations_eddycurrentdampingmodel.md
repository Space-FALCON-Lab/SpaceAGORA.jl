---
id: dynamics.perturbations_eddycurrentdampingmodel
label: EddyCurrentDampingModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: EddyCurrentDampingModel
  lines:
  - 1860
  - 1860
inputs:
- id: k_e
  type: Float64
  units: n/a
  required: true
  description: Field `k_e`.
- id: field_model
  type: Symbol
  units: n/a
  required: true
  description: Field `field_model`.
- id: igrf_year
  type: Float64
  units: n/a
  required: true
  description: Field `igrf_year`.
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
  type: EddyCurrentDampingModel
  units: n/a
  description: Constructed `EddyCurrentDampingModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# EddyCurrentDampingModel

## Purpose
Effector applying eddy-current damping torque from the spacecraft's rotation through Earth's magnetic field.

## Design & Implementation
Immutable with damping coefficient `k_e`, `field_model` of `:dipole` or `:igrf`, and `igrf_year`. The constructor requires positive finite `k_e`, validates the field model, requires a decimal year in `[1900, 2035)` for IGRF and warns once for years past 2030. Its `calcForceTorque` returns zero unless orientation is simulated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `k_e` | Float64 | n/a | yes | Field `k_e`. |
| in | `field_model` | Symbol | n/a | yes | Field `field_model`. |
| in | `igrf_year` | Float64 | n/a | yes | Field `igrf_year`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EddyCurrentDampingModel | n/a | — | Constructed `EddyCurrentDampingModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1875-1875`
- `callees` → [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1882-1882`
<!-- vulcan:connections:end -->

## Limitations
Earth-specific field models; no effect on translational state.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1860.
