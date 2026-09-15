---
id: core.reference_system_rtn_accel_to_inertial
label: rtn_accel_to_inertial
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtn_accel_to_inertial
  lines:
  - 698
  - 698
inputs:
- id: a_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `a_rtn`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `rtn_accel_to_inertial`.
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

# rtn_accel_to_inertial

## Purpose
Rotates an RTN acceleration command into the inertial frame so it can be applied as a force in the propagator.

## Design & Implementation
Multiplies the RTN DCM from the target state by the command vector. Because acceleration is a free vector, no transport or origin term applies; the DCM alone suffices.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a_rtn` | Any | n/a | yes | Positional argument `a_rtn`. |
| in | `r_target_ii` | Any | n/a | yes | Positional argument `r_target_ii`. |
| in | `v_target_ii` | Any | n/a | yes | Positional argument `v_target_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `rtn_accel_to_inertial`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_rtn_to_inertial_relative_state|rtn_to_inertial_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:694-694`
- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:18-18`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:18-18`

**Downstream**

- `callees` → [[core.reference_system_rtn_dcm_from_inertial|rtn_dcm_from_inertial]] · `callers` · call · `src/core/interfaces/reference_system.jl:699-699`
<!-- vulcan:connections:end -->

## Limitations
Recomputes the DCM on every call, including its two normalisations; a controller issuing commands at high rate could cache the DCM per tick.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 698.
