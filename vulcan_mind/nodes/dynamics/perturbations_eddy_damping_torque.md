---
id: dynamics.perturbations_eddy_damping_torque
label: eddy_damping_torque
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: eddy_damping_torque
  lines:
  - 1886
  - 1886
inputs:
- id: k_e
  type: Float64
  units: n/a
  required: true
  description: Positional argument `k_e`.
- id: B_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `B_body`.
- id: omega_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `ω_body`.
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
  description: Return value of `eddy_damping_torque`.
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

# eddy_damping_torque

## Purpose
The eddy-current damping law: the body-frame torque that opposes rotation of a conducting structure through the local magnetic field.

## Theory & Math
$$
\vec{\tau} = k_e\, \vec{B} \times (\vec{B} \times \vec{\omega})
$$

## Design & Implementation
Returns `k_e B × (B × ω)` with `k_e` in N·m·s/T², `B` in tesla and `ω` in rad/s, all in body coordinates. The double cross product projects the angular rate onto the plane perpendicular to the field, so rotation about the field direction is undamped. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `k_e` | Float64 | n/a | yes | Positional argument `k_e`. |
| in | `B_body` | SVector{3, Float64} | n/a | yes | Positional argument `B_body`. |
| in | `omega_body` | SVector{3, Float64} | n/a | yes | Positional argument `ω_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `eddy_damping_torque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_eddycurrentdampingmodel|EddyCurrentDampingModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1882-1882`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1917-1917`
- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1918-1918`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1910-1910`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1925-1925`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- `callees` → [[dynamics.perturbations__magnetic_field_inertial|_magnetic_field_inertial]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1911-1911`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- `callees` → [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1913-1913`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1925-1925`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1925-1925`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1907-1907`
<!-- vulcan:connections:end -->

## Limitations
A scalar `k_e` models an isotropic conductor; a real structure's eddy response is a tensor that depends on geometry and conductivity distribution.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1886.
