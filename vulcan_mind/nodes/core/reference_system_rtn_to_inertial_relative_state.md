---
id: core.reference_system_rtn_to_inertial_relative_state
label: rtn_to_inertial_relative_state
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtn_to_inertial_relative_state
  lines:
  - 677
  - 677
inputs:
- id: r_rel_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_rel_rtn`.
- id: v_rel_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_rel_rtn`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `rtn_to_inertial_relative_state`.
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

# rtn_to_inertial_relative_state

## Purpose
Reconstructs a chaser's inertial position and velocity from its RTN relative state and the target's inertial state, the exact inverse of the forward transform.

## Design & Implementation
Builds the DCM and rate, rotates the relative position by `C` and adds the target position, and rotates the relative velocity plus the transport term `n k̂ × r_rel` by `C` before adding the target velocity. Returns a tuple of two `SVector{3}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_rel_rtn` | Any | n/a | yes | Positional argument `r_rel_rtn`. |
| in | `v_rel_rtn` | Any | n/a | yes | Positional argument `v_rel_rtn`. |
| in | `r_target_ii` | Any | n/a | yes | Positional argument `r_target_ii`. |
| in | `v_target_ii` | Any | n/a | yes | Positional argument `v_target_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `rtn_to_inertial_relative_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:672-672`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system__rtn_rate_rad_s|_rtn_rate_rad_s]] · `callers` · call · `src/core/interfaces/reference_system.jl:684-684`
- `callees` → [[core.reference_system_rtn_accel_to_inertial|rtn_accel_to_inertial]] · `callers` · call · `src/core/interfaces/reference_system.jl:694-694`
- `callees` → [[core.reference_system_rtn_dcm_from_inertial|rtn_dcm_from_inertial]] · `callers` · call · `src/core/interfaces/reference_system.jl:683-683`
<!-- vulcan:connections:end -->

## Limitations
Shares the pure-cross-track-rotation assumption of the forward transform; round-tripping through both is exact only under that assumption.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 677.
