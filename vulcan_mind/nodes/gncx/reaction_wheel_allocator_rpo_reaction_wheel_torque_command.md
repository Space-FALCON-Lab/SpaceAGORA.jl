---
id: gncx.reaction_wheel_allocator_rpo_reaction_wheel_torque_command
label: rpo_reaction_wheel_torque_command
kind: function
source:
  file: src/gnc/control/rpo_mpc/reaction_wheel_allocator.jl
  symbol: rpo_reaction_wheel_torque_command
  lines:
  - 2
  - 15
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the RPO reaction-wheel allocator with
    the control model and the spacecraft attitude state.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rw_torque_body
  type: SVector{3,Float64}
  units: N m
  description: Saturated body-frame reaction-wheel torque command for attitude regulation
    during proximity operations.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# rpo_reaction_wheel_torque_command

## Purpose
`rpo_reaction_wheel_torque_command` produces the reaction-wheel torque that holds chaser attitude while the translational MPC law flies the approach. It is a proportional-derivative regulator on the attitude quaternion and body rate rather than a predictive law, because attitude regulation during proximity operations is a much faster and simpler loop than the relative translation problem.

## Model & Assumptions
The gains `attitude_kp` and `rate_kd` live on the `RPOMPCControlModel`. When both are exactly zero the routine returns a zero torque immediately, which is how attitude control is switched off without introducing a separate flag. The error signal is the vector part of the attitude quaternion, sign-corrected so the scalar part is non-negative; that correction selects the shorter of the two equivalent rotations and prevents the well-known unwinding behaviour where a regulator drives the body the long way round.

## Design & Implementation
The commanded torque is the negated sum of the proportional term on the sign-corrected quaternion vector and the derivative term on the body rate. Saturation is applied element-wise with `clamp.` against `max_rw_torque_nm`, and only when that limit is finite, so the default of infinity leaves the command untouched. The result is returned as a static three-vector, and `calcControlEffect!` in `rpo_mpc_control_model.jl` adds it to the thruster torque before caching the combined wrench in the `RPOHeldActuation` record.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the RPO reaction-wheel allocator with the control model and the spacecraft attitude state. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rw_torque_body` | SVector{3,Float64} | N m | — | Saturated body-frame reaction-wheel torque command for attitude regulation during proximity operations. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:27-27`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Element-wise clamping distorts the torque direction under saturation, unlike the norm-preserving scaling used by the momentum manager. Wheel momentum, wheel speed limits, and wheel geometry are not modelled here, so the command is a body torque request rather than a per-wheel allocation despite the file name. The regulator has no integral action and therefore leaves a steady-state error under a constant disturbance torque.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/reaction_wheel_allocator.jl:2-15`.
