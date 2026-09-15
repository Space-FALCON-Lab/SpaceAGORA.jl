---
id: gncx.rpo_control_types_rpompccontrolmodel
label: RPOMPCControlModel
kind: struct
source:
  file: src/gnc/control/rpo_mpc/rpo_control_types.jl
  symbol: RPOMPCControlModel
  lines:
  - 10
  - 21
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace that declares and exports the RPO model-predictive
    control effector type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rpo_model
  type: RPOMPCControlModel
  units: n/a
  description: Mutable control-effector record holding chaser and target indices,
    thruster geometry, controller, plan buffer, held actuation, and attitude gains.
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
# RPOMPCControlModel

## Purpose
`RPOMPCControlModel` is the state record for the rendezvous-and-proximity-operations control effector. It is a mutable keyword-constructed subtype of `AbstractControlEffectorModel`, so it participates in the four-function control contract exported by `ControlHooks` while carrying everything the MPC loop needs between updates.

## Model & Assumptions
The record separates three concerns. Identity is the chaser and target spacecraft indices, which let the effector pull both states out of the shared state vector and form the relative state. Actuation is the `SixAxisThrusterModel` geometry plus the attitude gains `attitude_kp` and `rate_kd` and the wheel saturation `max_rw_torque_nm`, whose defaults of zero and infinity mean attitude control is off unless configured. Timing and planning are the control step `control_dt_s` and the `RPOPlanBuffer`, which carries the guidance plan together with the time at which it was published so the reference preview can be indexed by elapsed plan time rather than absolute time.

## Design & Implementation
The `controller` field is deliberately typed `Any` and defaults to nothing, because the concrete `RpoLQMPCController` is built later by `init_rpo_lqmpc` and holds an OSQP model that cannot be sensibly default-constructed; `calcControlEffect!` returns immediately when it is still nothing. The companion `RPOHeldActuation` in the same file is the zero-order hold: it caches the inertial force, the body torque, the six thruster forces, and the reaction-wheel torque produced by the last control update, so the dynamics right-hand side can read a consistent command at every integrator stage between control ticks rather than re-solving the quadratic program inside the derivative evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace that declares and exports the RPO model-predictive control effector type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rpo_model` | RPOMPCControlModel | n/a | — | Mutable control-effector record holding chaser and target indices, thruster geometry, controller, plan buffer, held actuation, and attitude gains. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.rpo_control_types_rpoheldactuation|RPOHeldActuation]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:16-16`
- `callees` → [[gnc.rpo_plan_buffer_rpoplanbuffer|RPOPlanBuffer]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:15-15`
- `callees` → [[vehicle.thruster_models_sixaxisthrustermodel|SixAxisThrusterModel]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
Because the type is mutable and holds solver state, sharing one instance across spacecraft or across concurrent simulations is unsafe. The `Any`-typed controller field defeats specialisation at the call site. The zero-order hold means the applied command lags the true state by up to one control step, and nothing in the type enforces that `control_dt_s` is consistent with the discretisation step used to build the controller.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_control_types.jl:10-21`, with the held-actuation record at lines 2-7.
