---
id: dynamics.aerodynamic_wrench_models_aerodynamiccoefficientconstant
label: AerodynamicCoefficientConstant
kind: struct
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: AerodynamicCoefficientConstant
  lines:
  - 57
  - 57
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
  type: AerodynamicCoefficientConstant
  units: n/a
  description: Constructed `AerodynamicCoefficientConstant` (keyword constructor via
    @kwdef).
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

# AerodynamicCoefficientConstant

## Purpose
Marker effector type selecting the constant-coefficient drag model (`CD` linear in incidence from 0.8 to 2.2, zero lift) for spacecraft aerodynamics.

## Design & Implementation
An empty `@kwdef struct` subtyping `AbstractForceTorqueModel`. Dispatch on it selects `wrench`/`wrench_caching!` methods that call `_aero_pure_wrench(:constant, ...)`, `environment_requirements` returning `planet_frame=true, atmosphere=true`, and `solver_partition` returning `:implicit`. It has no `per_link_atmosphere` field, so per-link sampling for it can only be enabled through the process-wide `set_per_link_atmosphere!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerodynamicCoefficientConstant | n/a | — | Constructed `AerodynamicCoefficientConstant` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The legacy `calcForceTorque` method for this type references undefined variables (`vel_pp_rw_hat`, `q`, `args`, `rot_body_to_inertial`, `CS_body`, `drag_ii`, `drag_pp`) and will throw `UndefVarError` if invoked; only the `wrench` path is usable. Lacking `per_link_atmosphere` means it cannot opt in per instance.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 57.
