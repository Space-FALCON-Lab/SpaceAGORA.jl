---
id: analysis.scenario_builders__base_gravity_effector
label: _base_gravity_effector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _base_gravity_effector
  lines:
  - 20
  - 20
inputs:
- id: gravity_model
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `gravity_model`.
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
  type: Union{InverseSquaredGravityModel, InverseSquaredJ2GravityModel}
  units: n/a
  description: Return value of `_base_gravity_effector`. Returns `InverseSquaredGravityModel()`
    or `InverseSquaredJ2GravityModel()`.
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

# _base_gravity_effector

## Purpose
Maps a scenario's simple gravity choice onto the effector object the dynamics model uses when no spherical-harmonics file is configured.

## Design & Implementation
Returns `InverseSquaredGravityModel()` for `:inverse_squared` and `InverseSquaredJ2GravityModel()` for `:inverse_squared_j2`, raising `ArgumentError` otherwise. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gravity_model` | Symbol | n/a | yes | Positional argument `gravity_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{InverseSquaredGravityModel, InverseSquaredJ2GravityModel} | n/a | — | Return value of `_base_gravity_effector`. Returns `InverseSquaredGravityModel()` or `InverseSquaredJ2GravityModel()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:143-143`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[envana.env_gravity_models_inversesquaredj2gravitymodel|InverseSquaredJ2GravityModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:24-24`
- `callees` → [[environment.gravity_models_inversesquaredgravitymodel|InverseSquaredGravityModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:22-22`
<!-- vulcan:connections:end -->

## Limitations
It is bypassed entirely when `gravity_harmonics_degree` is positive, so a manifest that sets both a harmonics file and a gravity model symbol has the symbol silently ignored.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 20.
