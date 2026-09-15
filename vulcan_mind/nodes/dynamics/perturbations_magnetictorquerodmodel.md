---
id: dynamics.perturbations_magnetictorquerodmodel
label: MagneticTorqueRodModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: MagneticTorqueRodModel
  lines:
  - 1974
  - 1974
inputs:
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
  type: MagneticTorqueRodModel
  units: n/a
  description: Constructed `MagneticTorqueRodModel`.
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

# MagneticTorqueRodModel

## Purpose
Effector producing the torque from the spacecraft's magnets interacting with Earth's field.

## Design & Implementation
Immutable with `field_model` and `igrf_year`, validated as in the eddy model. Its torque is the total dipole moment crossed with the body-frame field, and zero unless orientation is simulated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `field_model` | Symbol | n/a | yes | Field `field_model`. |
| in | `igrf_year` | Float64 | n/a | yes | Field `igrf_year`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MagneticTorqueRodModel | n/a | — | Constructed `MagneticTorqueRodModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1991-1991`
- `callees` → [[dynamics.perturbations__magnetic_field_inertial|_magnetic_field_inertial]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1998-1998`
<!-- vulcan:connections:end -->

## Limitations
Earth-specific fields; the magnets' positions are ignored.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1974.
