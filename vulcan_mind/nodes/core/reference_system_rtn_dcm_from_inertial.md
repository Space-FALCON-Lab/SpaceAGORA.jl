---
id: core.reference_system_rtn_dcm_from_inertial
label: rtn_dcm_from_inertial
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtn_dcm_from_inertial
  lines:
  - 623
  - 623
inputs:
- id: r_target_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_target_ii`.
- id: v_target_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_target_ii`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `rtn_dcm_from_inertial`.
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

# rtn_dcm_from_inertial

## Purpose
Builds the rotation whose columns are the target's radial, along-track and cross-track unit vectors in inertial coordinates, the basis every HCW relative-motion calculation uses.

## Design & Implementation
Converts inputs to static vectors, raises `ArgumentError` if position or angular momentum is at or below machine epsilon, and forms `r̂`, `n̂ = h / |h|` and `t̂ = n̂ × r̂` normalised. Returns them concatenated as an `SMatrix{3,3}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_target_ii` | Any | n/a | yes | Positional argument `r_target_ii`. |
| in | `v_target_ii` | Any | n/a | yes | Positional argument `v_target_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `rtn_dcm_from_inertial`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:659-659`
- [[core.reference_system_rotate_vector_by_quaternion|rotate_vector_by_quaternion]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:617-617`
- [[core.reference_system_rtn_accel_to_inertial|rtn_accel_to_inertial]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:699-699`
- [[core.reference_system_rtn_to_inertial_relative_state|rtn_to_inertial_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:683-683`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing enforces that `v` is not parallel to `r` beyond the angular-momentum norm check; a nearly radial velocity gives a valid but numerically poor frame.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 623.
