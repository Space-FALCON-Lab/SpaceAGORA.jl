---
id: gnc.propulsive_maneuvers_calccontrolforcetorque
label: calcControlForceTorque
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: calcControlForceTorque
  lines:
  - 373
  - 373
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `calcControlForceTorque`.
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

# calcControlForceTorque

## Purpose
Evaluates the instantaneous control force and torque produced by a `BaseThrusterModel` for spacecraft `i` at time `t`, as the inertial-frame vectors consumed by the dynamics right-hand side.

## Design & Implementation
Returns a pair of zero `SVector{3, Float64}` when `i` is outside `1:length(controlModel.thrust)`. It obtains the window from `_effective_burn_window` and produces zero force outside the closed interval `[start_time, stop_time]`. Inside the window it takes the thrust magnitude from `_effective_thrust_isp`, builds `vel_vec` from `u.vel`, and computes `vel_mag = norm(vel_vec)`. Zero or non-finite speed, non-finite or non-positive thrust, or a non-finite direction each force a zero result. The thrust direction is `normalize(vel_vec)`, the prograde unit vector, multiplied by `+1.0` when `cos(direction_rad) >= 0.0` and `-1.0` otherwise, so zero radians is prograde and `π` is retrograde. Force is `thrust_mag * thrust_dir`; torque is returned as exactly zero, the source noting that offset and gimbaled thrusters are not yet modelled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `calcControlForceTorque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callees` → `callers` · feedback · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:386-386`

**Downstream**

- `callees` → [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- `callees` → [[gnc.propulsive_maneuvers__effective_burn_window|_effective_burn_window]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:378-378`
- `callees` → [[gnc.propulsive_maneuvers__effective_direction_rad|_effective_direction_rad]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:383-383`
- `callees` → [[gnc.propulsive_maneuvers__effective_thrust_isp|_effective_thrust_isp]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:380-380`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:400-400`
- `callees` → [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- `callees` → [[gnc.robot_arm_control_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:400-400`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:400-400`
- `callees` → [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- `callees` → [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- `callees` → [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
<!-- vulcan:connections:end -->

## Limitations
Thrust is confined to the velocity direction: only the sign of `cos(direction_rad)` survives, so any commanded out-of-plane or radial component is silently projected onto the along-track axis. The zero torque means thrust misalignment imparts no attitude disturbance at all, which understates the coupling a real offset thruster produces. The state is read as `u.vel`, whereas `calcControlEffect!` reads `u.sc[i].vel`, so the two entry points expect differently shaped state views. The window test uses no tolerance, unlike the 1e-9 second pad applied during scheduling.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 373.
