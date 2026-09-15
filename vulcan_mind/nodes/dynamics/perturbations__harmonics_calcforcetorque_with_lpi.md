---
id: dynamics.perturbations__harmonics_calcforcetorque_with_lpi
label: _harmonics_calcforcetorque_with_lpi
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_calcforcetorque_with_lpi
  lines:
  - 1509
  - 1509
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: AbstractVector{Float64}
  units: n/a
  required: true
  description: Positional argument `x`.
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: L_PI
  type: SMatrix{3, 3, Float64, 9}
  units: n/a
  required: true
  description: Positional argument `L_PI`.
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
  description: Return value of `_harmonics_calcforcetorque_with_lpi`.
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

# _harmonics_calcforcetorque_with_lpi

## Purpose
Evaluates the spherical-harmonics gravitational acceleration on one satellite using a supplied planet-frame rotation, the Pines-formulation core of the harmonics effector.

## Design & Implementation
Rotates the inertial position into the planet frame, obtains the satellite's scratch workspace, and forms the direction cosines `s, t, u` and the radial ratio. It builds the Helmholtz polynomial array `A` column by column using the precomputed `sqrt(2n+3)` factors, the real and imaginary sectoral recursions `R` and `I`, and accumulates the four Pines sums over the active orders per degree. The result is rotated back to inertial through the transpose and returned as force plus zero torque. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector{Float64} | n/a | yes | Positional argument `x`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `L_PI` | SMatrix{3, 3, Float64, 9} | n/a | yes | Positional argument `L_PI`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `_harmonics_calcforcetorque_with_lpi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1518-1518`
- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1663-1663`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- `callees` → [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1657-1657`
- `callees` → [[dynamics.perturbations__harmonics_workspace_for_sat_bang|_harmonics_workspace_for_sat!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1520-1520`
- `callees` → [[dynamics.perturbations__make_harmonics_scratch_workspace|_make_harmonics_scratch_workspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1672-1672`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1663-1663`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1663-1663`
<!-- vulcan:connections:end -->

## Limitations
Roughly 300 lines of recurrence with hand-indexed arrays; the `active_orders_by_degree` optimisation skips zero coefficients but makes the loop structure irregular, and no gravity-gradient torque is produced.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1509.
