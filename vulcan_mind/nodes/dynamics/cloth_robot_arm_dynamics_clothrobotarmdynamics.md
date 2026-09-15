---
id: dynamics.cloth_robot_arm_dynamics_clothrobotarmdynamics
label: ClothRobotArmDynamics
kind: module
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: ClothRobotArmDynamics
  lines:
  - 2
  - 2
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# ClothRobotArmDynamics

## Purpose
Module that realises a planned rigid robot arm (`RobotArmPlan`) as a compliant multibody chain of cloth-style bodies joined by spring-damper joints, so arm motion can be simulated standalone or coupled into a spacecraft state vector.

## Design & Implementation
Depends on `ClothMultibody` (compliant bodies, joints, actuators, steppers), `RobotArmPlanning` (plans, `robot_arm_plan_sample`, `cloth_fk`), and `Robotics`. Exports the reference-state sampler, the `ClothRobotArmSimulation` result type, model builders (`cloth_robot_arm_multibody`, `cloth_robot_arm_actuators`), state packing helpers (`cloth_robot_arm_initial_state`, `coupled_cloth_robot_arm_state_shape`, `initialize_coupled_cloth_robot_arm_state!`), the coupled RHS writer `assign_coupled_cloth_robot_arm_rhs!`, and `simulate_cloth_robot_arm_plan`. Quaternions use the scalar-last `(x, y, z, w)` layout with `_Q_IDENTITY = (0,0,0,1)`. `ClothRobotArmDynamicsMode` is an alias for `Symbol`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every link is modelled as a solid cylinder about its parent vector, ignoring any explicit inertia the link definition might carry. Default joint stiffness (5e3 N/m translation, 15 N m/rad rotation) and damping are hard-coded keyword defaults repeated in three functions. Rest quaternions are recomputed from forward kinematics at every RHS evaluation and every trajectory sample, which is expensive for long plans.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 2.
