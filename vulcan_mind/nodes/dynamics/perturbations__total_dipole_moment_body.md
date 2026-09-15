---
id: dynamics.perturbations__total_dipole_moment_body
label: _total_dipole_moment_body
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _total_dipole_moment_body
  lines:
  - 2023
  - 2023
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
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
  description: Return value of `_total_dipole_moment_body`.
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

# _total_dipole_moment_body

## Purpose
Sums the dipole moments of every magnet on every link into one body-frame moment for the torque-rod model.

## Design & Implementation
Nested loop over links and their `magnets`, accumulating `magnet.m` into a static vector. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_total_dipole_moment_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2062-2062`
- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2063-2063`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2055-2055`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2070-2070`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- `callees` → [[dynamics.perturbations__magnetic_field_inertial|_magnetic_field_inertial]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2056-2056`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- `callees` → [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2058-2058`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2070-2070`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- `callees` → [[dynx.coupled_perturbations_calculate_magnetic_torque|calculate_magnetic_torque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2064-2064`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2070-2070`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2052-2052`
<!-- vulcan:connections:end -->

## Limitations
Treats all magnets as acting at the centre of mass; their `location` fields are ignored, which is correct for a pure couple.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2023.
