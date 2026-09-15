---
id: gnc.robot_arm_control_calccontrolmassflowrate
label: calcControlMassFlowRate
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: calcControlMassFlowRate
  lines:
  - 211
  - 211
inputs:
- id: model
  type: RobotArmControlEffector
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Float64
  units: n/a
  description: Return value of `calcControlMassFlowRate`.
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

# calcControlMassFlowRate

## Purpose
Control-effector hook reporting propellant mass flow rate for the robot-arm controller. Since the arm is driven by electric joint actuators and consumes no propellant, the effector declares a zero flow so the mass state of the spacecraft is unaffected by arm motion.

## Design & Implementation
Signature `calcControlMassFlowRate(model::RobotArmControlEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64`, returning the literal `0.0` kg/s regardless of arguments. It overrides the generic `AbstractControlEffectorModel` fallback in `propulsive_maneuvers.jl`, which also returns `0.0`, purely to make the effector's behaviour explicit at its definition site. All five arguments are unused.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RobotArmControlEffector | n/a | yes | Positional argument `model`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `calcControlMassFlowRate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:400-400`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:395-395`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The method is functionally redundant with the abstract fallback and adds one more dispatch target to maintain. Any future arm actuator model that does consume a resource (for example cold-gas joint thrusters) would need this method rewritten rather than parameterised. The spacecraft index `i` is not checked, which is harmless only because the result is always zero.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 211.
