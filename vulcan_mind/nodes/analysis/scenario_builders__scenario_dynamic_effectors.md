---
id: analysis.scenario_builders__scenario_dynamic_effectors
label: _scenario_dynamic_effectors
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _scenario_dynamic_effectors
  lines:
  - 116
  - 116
inputs:
- id: cfg
  type: AbstractScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: cd_scale
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `cd_scale` (default `1.0`).
- id: cr_override
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Keyword argument `cr_override` (default `nothing`).
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
  type: Tuple
  units: n/a
  description: Return value of `_scenario_dynamic_effectors`. Returns `Tuple(effectors)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _scenario_dynamic_effectors

## Purpose
Assembles the tuple of force models a telemetry scenario flies with — gravity, optional N-body, optional SRP and optional scaled drag — from its manifest configuration.

## Design & Implementation
Gravity is either a `GravitationalHarmonicsModel` built from the manifest's degree, resolved order, file and the two environment-driven normalisation and J2-source choices, or the base effector; a positive degree without a file, or a missing file, raises `ArgumentError`. `NBodyGravityModel` is added when `nbody_bodies` is non-empty. SRP uses `srp_area_m2` if positive, else the bus reference area, and `cr_override` when given. Drag builds `AerodynamicCoefficientfM` with the configured incidence mode and wraps it in the scaled variant unless `cd_scale` is one within 1e-12.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | AbstractScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `cd_scale` | Float64 | n/a | no | Keyword argument `cd_scale` (default `1.0`). |
| in | `cr_override` | Union{Nothing, Float64} | n/a | no | Keyword argument `cr_override` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_scenario_dynamic_effectors`. Returns `Tuple(effectors)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:560-560`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:603-603`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:158-158`
- `callees` → [[analysis.scenario_builders__base_gravity_effector|_base_gravity_effector]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:143-143`
- `callees` → [[analysis.scenario_builders__harmonics_order|_harmonics_order]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:135-135`
- `callees` → [[analysis.scenario_builders__nbody_primary_name|_nbody_primary_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:151-151`
- `callees` → [[analysis.scenario_builders__telemetry_coefficients_normalized_for_scenario|_telemetry_coefficients_normalized_for_scenario]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:138-138`
- `callees` → [[analysis.scenario_builders__telemetry_j2_source_for_scenario|_telemetry_j2_source_for_scenario]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:139-139`
- `callees` → [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:170-170`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamiccoefficientfm|AerodynamicCoefficientfM]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:165-165`
- `callees` → [[dynamics.perturbations_gravitationalharmonicsmodel|GravitationalHarmonicsModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:133-133`
- `callees` → [[dynamics.perturbations_nbodygravitymodel|NBodyGravityModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:149-149`
- `callees` → [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:161-161`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:131-131`
<!-- vulcan:connections:end -->

## Limitations
Effectors are collected in a `Vector{Any}` and converted to a `Tuple`, so the returned tuple's type depends on the manifest and every distinct combination triggers a fresh specialisation of the dynamics right-hand side.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 116.
