---
id: dynamics.robot_arm_reaction_effector_calcforcetorque
label: calcForceTorque
kind: function
source:
  file: src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl
  symbol: calcForceTorque
  lines:
  - 27
  - 27
inputs:
- id: model
  type: RobotArmReactionEffector
  units: n/a
  required: true
  description: Positional argument `model`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  description: Return value of `calcForceTorque`. Returns `SVector{3, Float64}(0.0,
    0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)`.
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

# calcForceTorque

## Purpose

This `calcForceTorque(model::RobotArmReactionEffector, sc_view, p, sat_idx::Int)` method is the dispatch entry point through which the robot-arm reaction effector is polled by the coupled dynamics loop. It answers with the force and torque the arm applies to the spacecraft base, in the same tuple form every other force-torque effector in `DynamicEffectors` returns.

## Design & Implementation

The body is two guard clauses followed by a fixed return. First, `sat_idx == model.spacecraft_idx || return ...` short-circuits when the caller is integrating a satellite this effector does not own. Second, `model.plan === nothing && return ...` short-circuits when no `RobotArmPlan` has been attached. Both guards, and the fall-through, return the identical pair `SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)` — force in newtons and torque in newton-metres, both expressed in the body frame. The source comment marks the fall-through as a placeholder hook, noting that reaction wrench estimation properly belongs with the Cloth/RNEA multibody dynamics.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RobotArmReactionEffector | n/a | yes | Positional argument `model`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `calcForceTorque`. Returns `SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:265-265`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

No reaction wrench is ever computed: all three exit paths return zero, so attaching this effector changes nothing about the trajectory. The model's impedance gains, `joint_actuators` vector and `force_scale`/`torque_scale` are not consulted, and `robot_arm_plan_sample` is never invoked, so the plan itself is only used as a non-`nothing` presence check. Because the return type is fixed to `SVector{3, Float64}`, a future implementation cannot propagate dual numbers for forward-mode sensitivity without changing this signature.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl` line 27.
