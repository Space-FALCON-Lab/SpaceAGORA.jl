---
id: gnc.rpo_mpc_control_model__rpo_control_state_pos_vel
label: _rpo_control_state_pos_vel
kind: function
source:
  file: src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl
  symbol: _rpo_control_state_pos_vel
  lines:
  - 2
  - 2
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  type: SVector
  units: n/a
  description: Return value of `_rpo_control_state_pos_vel`. Returns `SVector{3, Float64}(u.sc[idx].pos),
    SVector{3, Float64}(u.sc[idx].vel)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _rpo_control_state_pos_vel

## Purpose
Pulls the inertial position and velocity of one spacecraft out of the integrator state for the rendezvous and proximity operations MPC controller, used for both the chaser and the target.

## Design & Implementation
Takes the state `u` and integer `idx` and returns the tuple `(SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel))`, in metres and metres per second respectively. Returning statically sized vectors keeps the downstream relative-state conversion, `inertial_to_rtn_relative_state`, and the acceleration frame rotations allocation-free inside the control callback.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_rpo_control_state_pos_vel`. Returns `SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:12-12`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It assumes the flat ComponentVector layout with named `pos` and `vel` fields, so it cannot be used under the partitioned gravity-backbone state and would throw on property access there. `idx` is unchecked, so an out-of-range chaser or target index surfaces as a bounds error from ComponentArrays. Conversion to `Float64` prevents use with dual numbers for differentiating the controller.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl` line 2.
