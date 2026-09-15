---
id: dynamics.aerodynamic_wrench_models_aerodynamiccoefficientnoballisticflight
label: AerodynamicCoefficientNoBallisticFlight
kind: struct
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: AerodynamicCoefficientNoBallisticFlight
  lines:
  - 106
  - 106
inputs:
- id: per_link_atmosphere
  type: Bool
  units: n/a
  required: false
  description: Field `per_link_atmosphere` (default `false`).
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
  type: AerodynamicCoefficientNoBallisticFlight
  units: n/a
  description: Constructed `AerodynamicCoefficientNoBallisticFlight` (keyword constructor
    via @kwdef).
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

# AerodynamicCoefficientNoBallisticFlight

## Purpose
Effector type intended for Newtonian-flow blunt-body aerodynamics (sphere-cone) rather than ballistic free-molecular flight.

## Design & Implementation
`@kwdef struct` with a single `per_link_atmosphere::Bool=false` field. Its `wrench`/`wrench_caching!` methods currently reuse `_aero_pure_wrench(:constant, ...)`, so at runtime it behaves identically to `AerodynamicCoefficientConstant`; the Newtonian coefficient function `aerodynamic_coefficient_no_ballistic_flight` is defined but never wired into the wrench path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `per_link_atmosphere` | Bool | n/a | no | Field `per_link_atmosphere` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerodynamicCoefficientNoBallisticFlight | n/a | — | Constructed `AerodynamicCoefficientNoBallisticFlight` (keyword constructor via @kwdef). |
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
The name promises Newtonian aerodynamics but the active code path is the constant-CD model. The legacy `calcForceTorque` method for this type contains a debugging `println` and references undefined variables (`vel_pp_rw_hat`, `q`, `args`, `CS_body`, `rot_body_to_inertial`), so it throws if called.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 106.
