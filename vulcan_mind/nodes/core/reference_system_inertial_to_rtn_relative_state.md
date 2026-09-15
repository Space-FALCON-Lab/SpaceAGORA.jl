---
id: core.reference_system_inertial_to_rtn_relative_state
label: inertial_to_rtn_relative_state
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: inertial_to_rtn_relative_state
  lines:
  - 653
  - 653
inputs:
- id: r_chaser_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_chaser_ii`.
- id: v_chaser_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_chaser_ii`.
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
  type: SVector{6,
  units: n/a
  description: Return value of `inertial_to_rtn_relative_state`.
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

# inertial_to_rtn_relative_state

## Purpose
Expresses a chaser's position and velocity relative to a target in the target's rotating RTN frame, the state the RPO controller and planner work in.

## Theory & Math
$$
r_{rtn} = C^\top (r_c - r_t),\qquad v_{rtn} = C^\top (v_c - v_t) - n\,\hat{k} \times r_{rtn}
$$

## Design & Implementation
Builds the RTN DCM and frame rate from the target state, subtracts the target from the chaser in inertial coordinates, rotates the relative position by the DCM transpose, and rotates the relative velocity likewise before subtracting `ω × r_rel` with `ω = n k̂` to account for the frame rotation. Returns a six-element `SVector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_chaser_ii` | Any | n/a | yes | Positional argument `r_chaser_ii`. |
| in | `v_chaser_ii` | Any | n/a | yes | Positional argument `v_chaser_ii`. |
| in | `r_target_ii` | Any | n/a | yes | Positional argument `r_target_ii`. |
| in | `v_target_ii` | Any | n/a | yes | Positional argument `v_target_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{6, | n/a | — | Return value of `inertial_to_rtn_relative_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__rtn_rate_rad_s|_rtn_rate_rad_s]] · `callees` → `callers` · feedback · `src/core/interfaces/reference_system.jl:648-648`
- [[gnc.rpo_guidance_hooks_build_rpo_plan|build_rpo_plan]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:11-11`
- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:91-91`
- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:14-14`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:14-14`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:11-11`

**Downstream**

- `callees` → [[core.reference_system__rtn_rate_rad_s|_rtn_rate_rad_s]] · `callers` · call · `src/core/interfaces/reference_system.jl:660-660`
- `callees` → [[core.reference_system_rtn_dcm_from_inertial|rtn_dcm_from_inertial]] · `callers` · call · `src/core/interfaces/reference_system.jl:659-659`
- `callees` → [[core.reference_system_rtn_to_inertial_relative_state|rtn_to_inertial_relative_state]] · `callers` · call · `src/core/interfaces/reference_system.jl:672-672`
<!-- vulcan:connections:end -->

## Limitations
The frame is assumed to rotate purely about its cross-track axis, which neglects the rate of change of the orbit plane; correct for two-body motion and an approximation under perturbations.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 653.
