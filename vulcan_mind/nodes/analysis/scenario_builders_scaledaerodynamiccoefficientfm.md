---
id: analysis.scenario_builders_scaledaerodynamiccoefficientfm
label: ScaledAerodynamicCoefficientfM
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: ScaledAerodynamicCoefficientfM
  lines:
  - 88
  - 88
inputs:
- id: model
  type: AerodynamicCoefficientfM
  units: n/a
  required: true
  description: Field `model`.
- id: cd_scale
  type: Float64
  units: n/a
  required: true
  description: Field `cd_scale`.
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
  type: ScaledAerodynamicCoefficientfM
  units: n/a
  description: Constructed `ScaledAerodynamicCoefficientfM`.
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

# ScaledAerodynamicCoefficientfM

## Purpose
Wraps the free-molecular aerodynamic model so a telemetry study can scale drag by a constant factor when fitting simulated decay to flight.

## Design & Implementation
An immutable subtype of `AbstractForceTorqueModel` holding the wrapped `AerodynamicCoefficientfM` and a positive `cd_scale`, enforced by an inner constructor that raises `ArgumentError` on non-positive values. `calcForceTorque` delegates to the wrapped model and multiplies both force and torque by the scale. It declares `environment_requirements` with atmosphere and planet frame so the engine's density-without-aero diagnostic does not misfire, and marks itself thread-safe.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerodynamicCoefficientfM | n/a | yes | Field `model`. |
| in | `cd_scale` | Float64 | n/a | yes | Field `cd_scale`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ScaledAerodynamicCoefficientfM | n/a | — | Constructed `ScaledAerodynamicCoefficientfM`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:170-170`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:104-104`
- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- `callees` → [[simulation.setup__dynamic_effector_threadsafe|_dynamic_effector_threadsafe]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:97-97`
<!-- vulcan:connections:end -->

## Limitations
Scaling the torque by the drag factor is a modelling shortcut; a true coefficient change would alter the pressure distribution and hence the centre of pressure, which this does not represent.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 88.
