---
id: gncx.momentum_manager_magneticmomentummanagermodel
label: MagneticMomentumManagerModel
kind: struct
source:
  file: src/gnc/control/momentum_manager.jl
  symbol: MagneticMomentumManagerModel
  lines:
  - 46
  - 60
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace that declares and exports the magnetic momentum
    manager as a control effector model.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: held_dipole
  type: SVector{3,Float64}
  units: A m^2
  description: Held magnetorquer dipole command, saturated at the rod capability,
    together with the resulting body torque.
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
# MagneticMomentumManagerModel

## Purpose
`MagneticMomentumManagerModel` is the control effector that unloads reaction-wheel momentum with magnetic torque rods. It is a mutable keyword-constructed model that holds both its configuration and its discrete internal state, so a single instance carries the wheel-momentum accumulator across the whole run.

## Theory & Math
The unloading law is $\mathbf{m} = \mu\,\dfrac{\mathbf{h}_w \times \mathbf{B}_b}{\|\mathbf{B}_b\|^2}$, giving a rod torque $\boldsymbol\tau = \mathbf{m} \times \mathbf{B}_b$ whose component along $\mathbf{h}_w$ is dissipative. The accumulator advances as $\mathbf{h}_w \leftarrow \mathbf{h}_w - \boldsymbol\tau_{\text{cmd}}\,\Delta t$.

## Model & Assumptions
Momentum management uses the standard cross-product control law: given the accumulated wheel momentum and the body-frame magnetic field, the commanded dipole is proportional to the cross product of the two divided by the squared field magnitude. Wheel momentum is not measured but integrated: each discrete update subtracts the commanded attitude torque times the elapsed step from the accumulator, seeded on the first tick from `h_wheels_0`. Both the commanded attitude torque and the inertial magnetic field are supplied as user callables, `commanded_torque` and `b_field_ii`, which keeps the manager independent of any particular attitude controller or field model.

## Design & Implementation
State updates happen only inside `calcControlEffect!`, which returns immediately unless the satellite index matches and both callables are present. The inertial field is rotated into the body frame through the attitude quaternion using `rot`. A degenerate field, where the squared magnitude is not strictly positive, zeroes both the dipole and the held torque instead of dividing by zero. The dipole is scaled down uniformly when its norm exceeds `m_max_am2`, which preserves the command direction under saturation. `calcControlForceTorque` then returns zero force and the held torque, and a `ticks` counter records how many discrete updates were applied. The file deliberately defines no `calcControlMassFlowRate` method, because torque rods consume no propellant and the abstract fallback in `propulsive_maneuvers.jl` already returns zero; a looser method here would be ambiguous against that fallback.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace that declares and exports the magnetic momentum manager as a control effector model. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `held_dipole` | SVector{3,Float64} | A m^2 | — | Held magnetorquer dipole command, saturated at the rod capability, together with the resulting body torque. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/momentum_manager.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The wheel momentum is open-loop integrated from commanded torque, so any mismatch between commanded and actual wheel torque accumulates without correction. The dipole is computed from the instantaneous field only, giving no guarantee of unloading along the field-parallel direction. Reading `h_wheels` or `held_dipole_am2` back after a run requires `run_simulation(args; isolate_state=false)`, because the default deep-copies the configuration and advances a copy.

## Provenance
Mapped from `src/gnc/control/momentum_manager.jl:46-60`, with the update law at lines 63-102.
