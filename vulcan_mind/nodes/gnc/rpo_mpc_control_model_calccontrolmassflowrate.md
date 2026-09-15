---
id: gnc.rpo_mpc_control_model_calccontrolmassflowrate
label: calcControlMassFlowRate
kind: function
source:
  file: src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl
  symbol: calcControlMassFlowRate
  lines:
  - 51
  - 51
inputs:
- id: model
  type: RPOMPCControlModel
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
Computes the propellant mass flow rate in kilograms per second implied by the currently held six-axis thruster command, so the mass state of the chaser depletes consistently with the control effort.

## Design & Implementation
Returns `0.0` for any index other than `model.chaser_idx`. For the chaser it loops `@inbounds` over the six thruster slots and, for each positive thrust `model.held.thruster_forces_n[j]`, subtracts `thrust / (model.thrusters.isp_s[j] * 9.80665)` from the accumulator. The result is therefore negative or zero, matching the sign convention that mass flow is added directly to the mass derivative. The constant 9.80665 is standard gravity in metres per second squared, converting specific impulse in seconds into effective exhaust velocity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOMPCControlModel | n/a | yes | Positional argument `model`. |
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
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:395-395`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Thrusts that are negative or exactly zero are skipped, so a solver that produced a negative entry contributes no propellant cost and that asymmetry can bias the optimiser. The loop count is hard-coded to six, so a thruster set of any other size is silently truncated or reads out of range under `@inbounds`. A zero or missing `isp_s` entry produces a division by zero and an infinite flow rate rather than an error.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl` line 51.
