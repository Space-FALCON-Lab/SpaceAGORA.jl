---
id: dynamics.perturbations__lvlh_cascade_torque
label: _lvlh_cascade_torque
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _lvlh_cascade_torque
  lines:
  - 2174
  - 2174
inputs:
- id: model
  type: LVLHCascadeAttitudeControlModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_ii`.
- id: q_ib
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Positional argument `q_ib`.
- id: w_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `w_body`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_lvlh_cascade_torque`.
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

# _lvlh_cascade_torque

## Purpose
Computes the body torque from a two-loop LVLH attitude controller: an outer attitude-error loop producing a rate command, and an inner rate loop producing torque.

## Design & Implementation
Builds the LVLH triad from position and velocity, forms the LVLH-to-body rotation from the current quaternion, and multiplies by the inverse of the commanded LVLH-to-body quaternion to get the error rotation. It extracts the error quaternion by the largest-trace branch, forms the outer-loop rate command from `k_out` times the error vector clamped to `w_max`, then the torque as `k_rate` times the rate error plus feedforward, saturated at `tau_max`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | LVLHCascadeAttitudeControlModel | n/a | yes | Positional argument `model`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Positional argument `vel_ii`. |
| in | `q_ib` | SVector{4, Float64} | n/a | yes | Positional argument `q_ib`. |
| in | `w_body` | SVector{3, Float64} | n/a | yes | Positional argument `w_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_lvlh_cascade_torque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_lvlhcascadeattitudecontrolmodel|LVLHCascadeAttitudeControlModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2168-2168`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2195-2195`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2263-2263`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2263-2263`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2263-2263`
<!-- vulcan:connections:end -->

## Limitations
Per-axis saturation of both loops means the commanded torque direction is distorted near limits; gains are constant with no gain scheduling on inertia.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2174.
